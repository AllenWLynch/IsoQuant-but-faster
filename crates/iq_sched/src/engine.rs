//! The coordinator: two-phase run with memory-budget admission.
//!
//! One coordinator (the calling thread, never a pool worker) owns all scheduling state and
//! a priority queue of ready tasks; tasks run on one rayon pool and report back over a
//! channel. Inputs move into tasks and results move back, so there is no shared mutable
//! state and no lock beyond the channel.
//!
//! Ordering of ready tasks:
//! - class "now" (barrier, merge, spill): short, memory-neutral; admitted ahead of everything
//!   and without a memory check (only the running cap);
//! - class "pass A": largest chromosome first. Pass A is cheap relative to collect and
//!   unlocks work (and every chromosome must pass the barrier before Phase M), so it goes
//!   ahead of collect; pass A of one chromosome overlaps collect of the chromosomes already
//!   done (vtg_proto finding #11);
//! - class "work" (collect, construct): largest estimated cost first (LPT), strict
//!   head-of-line admission against the byte budget.
//!
//! At most `threads` tasks run at once; a task may use more workers through nested rayon
//! parallelism (a giant block's per-read work), so idle workers are never stranded.

use crate::admission::Admission;
use crate::pipeline::{Pipeline, ReadTotals};
use crate::profile::{cpu_seconds, Kind, Report, RssSampler, TaskRec};
use crate::units::{pack_target, schedule_units, SplitMode};
use iq_types::Block;
use std::collections::BinaryHeap;
use std::ops::Range;
use std::sync::{Arc, OnceLock};
use std::time::{Duration, Instant};

/// Scheduler options.
#[derive(Clone, Debug)]
pub struct EngineOpts {
    pub threads: usize,
    /// Total byte budget (before `Pipeline::fixed_bytes` is subtracted).
    pub budget: u64,
    pub split: SplitMode,
    /// Sample RSS every 100 ms (cheap; always on for `--profile`).
    pub sample_rss: bool,
    /// Log scheduling decisions to stderr.
    pub verbose: bool,
}

impl Default for EngineOpts {
    fn default() -> Self {
        EngineOpts { threads: 1, budget: u64::MAX / 4, split: SplitMode::Packed, sample_rss: false, verbose: false }
    }
}

/// A scheduling unit: consecutive blocks of one chromosome.
#[derive(Clone, Debug)]
pub struct UnitSpec {
    pub chr: usize,
    pub blocks: Range<usize>,
    pub reads: u64,
    pub cost: u64,
    pub est_collect: u64,
    pub est_construct: u64,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Task {
    PassA(usize),
    Collect(usize),
    Spill(usize),
    Barrier(usize),
    Construct(usize),
    Merge(usize),
}

#[derive(PartialEq, Eq)]
struct QItem {
    /// 2 = now (barrier/merge/spill), 1 = pass A, 0 = collect/construct.
    class: u8,
    cost: u64,
    seq: std::cmp::Reverse<u64>,
    task: TaskKey,
}

/// `Task` with a total order for the heap (ties never reach it: `seq` is unique).
#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
struct TaskKey(u8, usize);

impl TaskKey {
    fn of(t: Task) -> Self {
        match t {
            Task::PassA(i) => TaskKey(0, i),
            Task::Collect(i) => TaskKey(1, i),
            Task::Spill(i) => TaskKey(2, i),
            Task::Barrier(i) => TaskKey(3, i),
            Task::Construct(i) => TaskKey(4, i),
            Task::Merge(i) => TaskKey(5, i),
        }
    }
    fn task(self) -> Task {
        match self.0 {
            0 => Task::PassA(self.1),
            1 => Task::Collect(self.1),
            2 => Task::Spill(self.1),
            3 => Task::Barrier(self.1),
            4 => Task::Construct(self.1),
            _ => Task::Merge(self.1),
        }
    }
}

impl Ord for QItem {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        (self.class, self.cost, self.seq).cmp(&(o.class, o.cost, o.seq))
    }
}
impl PartialOrd for QItem {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}

enum Payload<P: Pipeline> {
    PassA(Result<Vec<Block>, String>),
    Collected(Result<Vec<P::Store>, String>, u64),
    Spilled(Result<Vec<P::Store>, String>, u64, u64),
    Barrier(Result<(Vec<P::Store>, ReadTotals), String>),
    /// A store was consumed inside a running construct task (not a task completion).
    Freed(u64),
    Constructed(Result<Vec<P::BlockOut>, String>),
    Merged(Result<P::ChrOut, String>),
}

struct Ev<P: Pipeline> {
    task: Option<Task>,
    start: Duration,
    end: Duration,
    payload: Payload<P>,
}

struct ChrState<P: Pipeline> {
    blocks: Arc<Vec<Block>>,
    units: Vec<usize>,
    stores: Vec<Option<P::Store>>,
    collected: usize,
    pending_spill: usize,
    barrier_queued: bool,
    outs: Vec<Option<P::BlockOut>>,
    constructed: usize,
    merged: Option<P::ChrOut>,
}

/// Runs `p` through the two-phase schedule. Returns the run report.
pub fn execute<P: Pipeline>(p: &P, opts: &EngineOpts) -> Result<Report, String> {
    let threads = opts.threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .thread_name(|i| format!("iqrs-w{i}"))
        .build()
        .map_err(|e| e.to_string())?;
    let sampler = opts.sample_rss.then(|| RssSampler::start(Duration::from_millis(100)));
    let baseline_rss = sampler.as_ref().map_or(0, |s| s.baseline);
    let t0 = Instant::now();
    let cpu0 = cpu_seconds();
    let budget = opts.budget.saturating_sub(p.fixed_bytes());

    let n_chr = p.n_chromosomes();
    let global: OnceLock<P::Global> = OnceLock::new();
    let mut co = Coord::<P> {
        p,
        opts,
        t0,
        adm: Admission::new(budget),
        cap: threads,
        queue: BinaryHeap::new(),
        seq: 0,
        chrs: (0..n_chr).map(|_| None).collect(),
        units: Vec::new(),
        unit_stores: Vec::new(),
        passa_done: 0,
        barriers_done: 0,
        merged: 0,
        totals: ReadTotals::default(),
        tasks: Vec::new(),
        error: None,
        phase_m: false,
    };
    for c in 0..n_chr {
        co.push(Task::PassA(c), p.chr_weight(c));
    }

    let (tx, rx) = crossbeam_channel::unbounded::<Ev<P>>();
    pool.in_place_scope(|scope| loop {
        if co.done() {
            break;
        }
        co.admit(scope, &tx, &global);
        if co.adm.running == 0 {
            if co.error.is_some() || co.done() {
                break;
            }
            if co.queue.is_empty() {
                co.error = Some("scheduler stalled: nothing running and nothing ready".into());
                break;
            }
            continue;
        }
        let ev = rx.recv().expect("channel closed while tasks run");
        co.on_event(ev, &global);
    });
    // Drain events of tasks that finished after the loop (none expected: running == 0).
    drop(tx);

    if let Some(e) = co.error.take() {
        return Err(e);
    }
    let outs: Vec<P::ChrOut> =
        co.chrs.iter_mut().map(|c| c.as_mut().and_then(|c| c.merged.take()).expect("merged")).collect();
    let fs = t0.elapsed();
    let totals = co.totals;
    pool.install(|| p.finish(outs, totals))?;
    co.tasks.push(TaskRec {
        kind: Kind::Finish,
        id: 0,
        start: fs,
        end: t0.elapsed(),
        est_bytes: 0,
        actual_bytes: 0,
        reads: 0,
    });

    let (peak_rss, rss_series) = sampler.map_or((0, vec![]), |s| s.stop());
    Ok(Report {
        threads,
        budget,
        wall: t0.elapsed(),
        cpu_seconds: cpu_seconds() - cpu0,
        peak_rss,
        baseline_rss,
        peak_accounted: co.adm.peak_accounted,
        over_budget_admissions: co.adm.over_budget_admissions,
        n_units: co.units.len(),
        tasks: co.tasks,
        rss_series,
    })
}

struct Coord<'a, P: Pipeline> {
    p: &'a P,
    opts: &'a EngineOpts,
    t0: Instant,
    adm: Admission,
    cap: usize,
    queue: BinaryHeap<QItem>,
    seq: u64,
    chrs: Vec<Option<ChrState<P>>>,
    units: Vec<UnitSpec>,
    /// Stores of a unit in flight between tasks (collect → spill → slots).
    unit_stores: Vec<Option<Vec<P::Store>>>,
    passa_done: usize,
    barriers_done: usize,
    merged: usize,
    totals: ReadTotals,
    tasks: Vec<TaskRec>,
    error: Option<String>,
    phase_m: bool,
}

impl<'a, P: Pipeline> Coord<'a, P> {
    fn done(&self) -> bool {
        self.merged == self.chrs.len()
    }

    fn push(&mut self, t: Task, cost: u64) {
        self.seq += 1;
        let class = match t {
            Task::Barrier(_) | Task::Merge(_) | Task::Spill(_) => 2,
            Task::PassA(_) => 1,
            Task::Collect(_) | Task::Construct(_) => 0,
        };
        self.queue.push(QItem { class, cost, seq: std::cmp::Reverse(self.seq), task: TaskKey::of(t) });
    }

    fn est(&self, t: Task) -> u64 {
        match t {
            Task::PassA(c) => self.p.est_pass_a(c),
            Task::Collect(u) => self.units[u].est_collect,
            Task::Construct(u) => self.units[u].est_construct,
            _ => 0,
        }
    }

    fn log(&self, msg: impl FnOnce() -> String) {
        if self.opts.verbose {
            eprintln!("[sched {:>8.3}s] {}", self.t0.elapsed().as_secs_f64(), msg());
        }
    }

    /// Admits ready tasks while they fit (strict head-of-line within each class).
    fn admit<'s>(
        &mut self,
        scope: &rayon::Scope<'s>,
        tx: &crossbeam_channel::Sender<Ev<P>>,
        global: &'s OnceLock<P::Global>,
    ) where
        'a: 's,
    {
        if self.error.is_some() {
            return;
        }
        while let Some(head) = self.queue.peek() {
            if self.adm.running >= self.cap {
                break;
            }
            let task = head.task.task();
            let est = self.est(task);
            if head.class == 2 {
                // memory-neutral, short: running cap only
                self.adm.running += 1;
            } else if !self.adm.try_admit(est) {
                break;
            }
            self.queue.pop();
            self.spawn(task, est, scope, tx, global);
        }
    }

    fn spawn<'s>(
        &mut self,
        task: Task,
        est: u64,
        scope: &rayon::Scope<'s>,
        tx: &crossbeam_channel::Sender<Ev<P>>,
        global: &'s OnceLock<P::Global>,
    ) where
        'a: 's,
    {
        let p: &'s P = self.p;
        let t0 = self.t0;
        let tx = tx.clone();
        self.log(|| format!("start {task:?} est {:.1} MiB; {:?}", est as f64 / 1048576.0, self.adm));
        let send = move |tx: &crossbeam_channel::Sender<Ev<P>>, start: Duration, payload: Payload<P>| {
            let _ = tx.send(Ev { task: Some(task), start, end: t0.elapsed(), payload });
        };
        match task {
            Task::PassA(c) => scope.spawn(move |_| {
                let st = t0.elapsed();
                send(&tx, st, Payload::PassA(p.pass_a(c)));
            }),
            Task::Collect(u) => {
                let spec = self.units[u].clone();
                let blocks = self.chrs[spec.chr].as_ref().unwrap().blocks.clone();
                scope.spawn(move |_| {
                    let st = t0.elapsed();
                    let mut stores = Vec::with_capacity(spec.blocks.len());
                    let mut res = Ok(());
                    for bi in spec.blocks.clone() {
                        match p.collect_block(spec.chr, bi, &blocks[bi]) {
                            Ok(s) => stores.push(s),
                            Err(e) => {
                                res = Err(e);
                                break;
                            }
                        }
                    }
                    let bytes = stores.iter().map(|s| p.store_bytes(s)).sum();
                    send(&tx, st, Payload::Collected(res.map(|_| stores), bytes));
                })
            }
            Task::Spill(u) => {
                let spec = self.units[u].clone();
                let mut stores = self.unit_stores[u].take().expect("stores to spill");
                scope.spawn(move |_| {
                    let st = t0.elapsed();
                    let before: u64 = stores.iter().map(|s| p.store_bytes(s)).sum();
                    let mut after = 0;
                    let mut res = Ok(());
                    for (s, bi) in stores.iter_mut().zip(spec.blocks.clone()) {
                        match p.spill(spec.chr, bi, s) {
                            Ok(b) => after += b,
                            Err(e) => {
                                res = Err(e);
                                break;
                            }
                        }
                    }
                    send(&tx, st, Payload::Spilled(res.map(|_| stores), before, after));
                })
            }
            Task::Barrier(c) => {
                let cs = self.chrs[c].as_mut().unwrap();
                let blocks = cs.blocks.clone();
                let mut stores: Vec<P::Store> = cs.stores.iter_mut().map(|s| s.take().expect("store")).collect();
                scope.spawn(move |_| {
                    let st = t0.elapsed();
                    let r = p.barrier_chr(c, &blocks, &mut stores).map(|t| (stores, t));
                    send(&tx, st, Payload::Barrier(r));
                })
            }
            Task::Construct(u) => {
                let spec = self.units[u].clone();
                let cs = self.chrs[spec.chr].as_mut().unwrap();
                let blocks = cs.blocks.clone();
                let stores: Vec<P::Store> =
                    spec.blocks.clone().map(|bi| cs.stores[bi].take().expect("store")).collect();
                scope.spawn(move |_| {
                    let st = t0.elapsed();
                    let g = global.get().expect("global decided before phase M");
                    let mut state = P::UnitState::default();
                    let mut outs = Vec::with_capacity(stores.len());
                    let mut res = Ok(());
                    for (store, bi) in stores.into_iter().zip(spec.blocks.clone()) {
                        let b = p.store_bytes(&store);
                        let r = p.construct_block(spec.chr, bi, &blocks[bi], store, &mut state, g);
                        let _ = tx.send(Ev { task: None, start: st, end: st, payload: Payload::Freed(b) });
                        match r {
                            Ok(o) => outs.push(o),
                            Err(e) => {
                                res = Err(e);
                                break;
                            }
                        }
                    }
                    send(&tx, st, Payload::Constructed(res.map(|_| outs)));
                })
            }
            Task::Merge(c) => {
                let cs = self.chrs[c].as_mut().unwrap();
                let blocks = cs.blocks.clone();
                let outs: Vec<P::BlockOut> = cs.outs.iter_mut().map(|o| o.take().expect("block out")).collect();
                scope.spawn(move |_| {
                    let st = t0.elapsed();
                    send(&tx, st, Payload::Merged(p.merge_chr(c, &blocks, outs)));
                })
            }
        }
    }

    fn fail(&mut self, e: String) {
        self.log(|| format!("error: {e}"));
        if self.error.is_none() {
            self.error = Some(e);
        }
    }

    fn on_event(&mut self, ev: Ev<P>, global: &OnceLock<P::Global>) {
        let task = ev.task;
        let est = task.map_or(0, |t| self.est(t));
        let is_now = matches!(task, Some(Task::Spill(_) | Task::Barrier(_) | Task::Merge(_)));
        // release the task's admission
        match &ev.payload {
            Payload::Freed(b) => {
                self.adm.release(*b);
                return;
            }
            Payload::Collected(_, bytes) => self.adm.finish(est, *bytes),
            _ if is_now => self.adm.running -= 1,
            _ => self.adm.finish(est, 0),
        }
        let task = task.expect("completion");
        match (task, ev.payload) {
            (Task::PassA(c), Payload::PassA(r)) => {
                let st = (ev.start, ev.end);
                match r {
                    Ok(blocks) => {
                        let reads: u64 = blocks.iter().map(|b| b.owned_reads).sum();
                        self.tasks.push(TaskRec {
                            kind: Kind::PassA,
                            id: c,
                            start: st.0,
                            end: st.1,
                            est_bytes: est,
                            actual_bytes: 0,
                            reads,
                        });
                        self.on_pass_a(c, blocks);
                    }
                    Err(e) => self.fail(format!("pass A {}: {e}", self.p.chr_name(c))),
                }
            }
            (Task::Collect(u), Payload::Collected(r, bytes)) => {
                let (reads, chr) = (self.units[u].reads, self.units[u].chr);
                self.tasks.push(TaskRec {
                    kind: Kind::Collect,
                    id: u,
                    start: ev.start,
                    end: ev.end,
                    est_bytes: est,
                    actual_bytes: bytes,
                    reads,
                });
                match r {
                    Ok(stores) => {
                        self.log(|| format!("collected unit {u} ({reads} reads): est {est} B, actual {bytes} B"));
                        if self.p.can_spill() && self.adm.resident > self.adm.budget / 2 {
                            self.unit_stores[u] = Some(stores);
                            self.chrs[chr].as_mut().unwrap().pending_spill += 1;
                            self.push(Task::Spill(u), 0);
                        } else {
                            self.place_stores(u, stores);
                        }
                        self.chrs[chr].as_mut().unwrap().collected += 1;
                        self.maybe_barrier(chr);
                    }
                    Err(e) => self.fail(format!("collect unit {u} ({}): {e}", self.p.chr_name(chr))),
                }
            }
            (Task::Spill(u), Payload::Spilled(r, before, after)) => {
                let chr = self.units[u].chr;
                match r {
                    Ok(stores) => {
                        self.tasks.push(TaskRec {
                            kind: Kind::Spill,
                            id: u,
                            start: ev.start,
                            end: ev.end,
                            est_bytes: before,
                            actual_bytes: after,
                            reads: 0,
                        });
                        self.adm.release(before);
                        self.adm.add_resident(after);
                        self.place_stores(u, stores);
                        self.chrs[chr].as_mut().unwrap().pending_spill -= 1;
                        self.maybe_barrier(chr);
                    }
                    Err(e) => self.fail(format!("spill unit {u}: {e}")),
                }
            }
            (Task::Barrier(c), Payload::Barrier(r)) => match r {
                Ok((stores, t)) => {
                    self.tasks.push(TaskRec {
                        kind: Kind::Barrier,
                        id: c,
                        start: ev.start,
                        end: ev.end,
                        est_bytes: 0,
                        actual_bytes: 0,
                        reads: t.total,
                    });
                    let cs = self.chrs[c].as_mut().unwrap();
                    for (slot, s) in cs.stores.iter_mut().zip(stores) {
                        *slot = Some(s);
                    }
                    self.totals += t;
                    self.barriers_done += 1;
                    self.maybe_phase_m(global);
                }
                Err(e) => self.fail(format!("barrier {}: {e}", self.p.chr_name(c))),
            },
            (Task::Construct(u), Payload::Constructed(r)) => {
                let spec = self.units[u].clone();
                self.tasks.push(TaskRec {
                    kind: Kind::Construct,
                    id: u,
                    start: ev.start,
                    end: ev.end,
                    est_bytes: est,
                    actual_bytes: 0,
                    reads: spec.reads,
                });
                match r {
                    Ok(outs) => {
                        let cs = self.chrs[spec.chr].as_mut().unwrap();
                        for (bi, o) in spec.blocks.clone().zip(outs) {
                            cs.outs[bi] = Some(o);
                        }
                        cs.constructed += 1;
                        if cs.constructed == cs.units.len() {
                            self.push(Task::Merge(spec.chr), 0);
                        }
                    }
                    Err(e) => self.fail(format!("construct unit {u} ({}): {e}", self.p.chr_name(spec.chr))),
                }
            }
            (Task::Merge(c), Payload::Merged(r)) => {
                self.tasks.push(TaskRec {
                    kind: Kind::Merge,
                    id: c,
                    start: ev.start,
                    end: ev.end,
                    est_bytes: 0,
                    actual_bytes: 0,
                    reads: 0,
                });
                match r {
                    Ok(o) => {
                        self.chrs[c].as_mut().unwrap().merged = Some(o);
                        self.merged += 1;
                    }
                    Err(e) => self.fail(format!("merge {}: {e}", self.p.chr_name(c))),
                }
            }
            _ => unreachable!("event/task mismatch"),
        }
    }

    fn place_stores(&mut self, u: usize, stores: Vec<P::Store>) {
        let spec = &self.units[u];
        let cs = self.chrs[spec.chr].as_mut().unwrap();
        for (bi, s) in spec.blocks.clone().zip(stores) {
            cs.stores[bi] = Some(s);
        }
    }

    fn on_pass_a(&mut self, c: usize, blocks: Vec<Block>) {
        let chr_reads: u64 = blocks.iter().map(|b| b.owned_reads).sum();
        let ranges = schedule_units(&blocks, self.opts.split, pack_target(chr_reads, self.opts.threads));
        let mut ids = Vec::with_capacity(ranges.len());
        for r in ranges {
            let bs = &blocks[r.clone()];
            let spec = UnitSpec {
                chr: c,
                reads: bs.iter().map(|b| b.owned_reads).sum(),
                cost: bs.iter().map(|b| self.p.est_cost(b)).sum(),
                est_collect: bs.iter().map(|b| self.p.est_collect(b)).sum(),
                est_construct: bs.iter().map(|b| self.p.est_construct(b)).max().unwrap_or(0),
                blocks: r,
            };
            let u = self.units.len();
            ids.push(u);
            self.push(Task::Collect(u), spec.cost);
            self.units.push(spec);
            self.unit_stores.push(None);
        }
        let n = blocks.len();
        self.log(|| format!("pass A {} done: {n} blocks, {} units", self.p.chr_name(c), ids.len()));
        self.chrs[c] = Some(ChrState {
            blocks: Arc::new(blocks),
            units: ids,
            stores: (0..n).map(|_| None).collect(),
            collected: 0,
            pending_spill: 0,
            barrier_queued: false,
            outs: (0..n).map(|_| None).collect(),
            constructed: 0,
            merged: None,
        });
        self.passa_done += 1;
        self.maybe_barrier(c);
    }

    fn maybe_barrier(&mut self, c: usize) {
        let cs = self.chrs[c].as_mut().unwrap();
        if !cs.barrier_queued && cs.collected == cs.units.len() && cs.pending_spill == 0 {
            cs.barrier_queued = true;
            self.push(Task::Barrier(c), 0);
        }
    }

    fn maybe_phase_m(&mut self, global: &OnceLock<P::Global>) {
        if self.phase_m || self.barriers_done < self.chrs.len() {
            return;
        }
        self.phase_m = true;
        let g = self.p.decide(self.totals);
        let _ = global.set(g);
        self.log(|| format!("barrier passed: totals {:?}; phase M starts", self.totals));
        for u in 0..self.units.len() {
            let cost = self.units[u].cost;
            self.push(Task::Construct(u), cost);
        }
        // chromosomes without units merge immediately
        for c in 0..self.chrs.len() {
            if self.chrs[c].as_ref().unwrap().units.is_empty() {
                self.push(Task::Merge(c), 0);
            }
        }
    }
}
