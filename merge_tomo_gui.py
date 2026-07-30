#!/usr/bin/env python3
"""
merge_tomo_gui.py -- desktop front end for merge_tomo.py.

Drop (or browse for) the west and east .xyz files, set the smoothing length,
press MERGE.  All numerics live in merge_tomo.py; this file is only the window,
so the GUI and the command line can never drift apart.

    python merge_tomo_gui.py

Drag-and-drop needs one optional package:

    pip install tkinterdnd2

Without it the app still works - the drop panels become click-to-browse.
"""

import os
import queue
import sys
import threading
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import merge_tomo as mt

# Optional drag-and-drop.
try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAVE_DND = True
except Exception:
    DND_FILES, TkinterDnD = None, None
    HAVE_DND = False

# Optional embedded preview.  The figure is built with Figure(), never pyplot,
# so merge_tomo.qc_plot's matplotlib.use("Agg") cannot disturb this canvas.
try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

BG      = "#f4f4f6"
PANEL   = "#ffffff"
EDGE    = "#c9ccd4"
HOT     = "#dbeafe"
ACCENT  = "#1f4e8c"
OK_C    = "#1a7f37"
WARN_C  = "#b45309"
BAD_C   = "#b42318"


# ============================================================================
# Pure logic - no Tk in here, so it can be unit-tested headlessly
# ============================================================================
def parse_drop(data):
    """
    Turn a tkinterdnd2 drop payload into a list of paths.

    Tk hands over a space-separated list in which any path containing spaces is
    wrapped in braces:  '{C:/My Data/west.xyz} C:/east.xyz'
    """
    out, buf, depth = [], "", 0
    for ch in data:
        if ch == "{":
            depth += 1
            if depth == 1:
                continue
        if ch == "}":
            depth -= 1
            if depth == 0:
                out.append(buf); buf = ""
                continue
        if ch == " " and depth == 0:
            if buf:
                out.append(buf); buf = ""
            continue
        buf += ch
    if buf:
        out.append(buf)
    return [p for p in out if p]


def summarise(lon, lat, vel):
    return (f"{len(lon)} nodes\n"
            f"lon {lon.min():.2f} .. {lon.max():.2f}\n"
            f"lat {lat.min():.2f} .. {lat.max():.2f}\n"
            f"vel {vel.min():.3f} .. {vel.max():.3f}")


def order_by_longitude(loaded):
    """Given [(path, data), ...] return them sorted west-most first."""
    return sorted(loaded, key=lambda pd: float(np.mean(pd[1][0])))


def default_centres(west, east):
    return (round((west[0].min() + west[0].max()) / 2, 2),
            round((east[0].min() + east[0].max()) / 2, 2))


def overlap_span(west, east):
    lo = max(west[0].min(), east[0].min())
    hi = min(west[0].max(), east[0].max())
    return lo, hi


def edge_quality(L, cW, cE, ov_lo, ov_hi):
    """
    Weights at both overlap edges plus a verdict.

    The blend has to be ~pure west where the overlap ends in the west and
    ~pure east where it ends in the east; otherwise the merged map steps at
    those longitudes.  Returns (w_at_west_edge, east_weight_at_east_edge,
    verdict) with verdict in {'ok', 'loose', 'none'}.
    """
    if ov_hi <= ov_lo or cE == cW:
        return None, None, "none"
    wlo = float(mt.west_weight(ov_lo, L, cW, cE))
    whi = float(mt.west_weight(ov_hi, L, cW, cE))
    worst = max(1.0 - wlo, whi)
    return wlo, 1.0 - whi, ("ok" if worst <= mt.EDGE_TOL + 1e-6 else "loose")


def slider_range(cW, cE, ov_lo, ov_hi):
    """Sensible slider bounds: comfortably either side of the recommended L."""
    Ls = mt.suggest_L(cW, cE, ov_lo, ov_hi)
    hi = max(2.5 * Ls, 1.0) if Ls else max(abs(cE - cW), 10.0)
    return 0.05, round(hi, 2)


def default_output(path_west):
    d = os.path.dirname(os.path.abspath(path_west)) if path_west else os.getcwd()
    return os.path.join(d, "merged.xyz")


def open_externally(path):
    """Best-effort 'show me this file' across platforms."""
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)                                # noqa: S606
        elif sys.platform == "darwin":
            import subprocess; subprocess.Popen(["open", path])
        else:
            import subprocess; subprocess.Popen(["xdg-open", path])
        return True
    except Exception:
        return False


# ============================================================================
# Drop / browse panel
# ============================================================================
class DropPanel(tk.Frame):
    def __init__(self, master, title, on_files):
        super().__init__(master, bg=PANEL, highlightbackground=EDGE,
                         highlightthickness=1, bd=0)
        self.on_files = on_files
        self.title = title

        self.lbl_title = tk.Label(self, text=title, bg=PANEL, fg=ACCENT,
                                  font=("Segoe UI", 11, "bold"))
        self.lbl_title.pack(anchor="w", padx=10, pady=(8, 0))

        hint = ("drop a .xyz file here" if HAVE_DND else "click to choose a .xyz file")
        self.lbl_path = tk.Label(self, text=hint, bg=PANEL, fg="#6b7280",
                                 font=("Segoe UI", 9), anchor="w",
                                 wraplength=300, justify="left")
        self.lbl_path.pack(anchor="w", padx=10, pady=(2, 0))

        self.lbl_info = tk.Label(self, text="", bg=PANEL, fg="#374151",
                                 font=("Consolas", 9), anchor="w",
                                 justify="left")
        self.lbl_info.pack(anchor="w", padx=10, pady=(6, 10))

        for w in (self, self.lbl_title, self.lbl_path, self.lbl_info):
            w.bind("<Button-1>", lambda _e: self.browse())
        if HAVE_DND:
            self.drop_target_register(DND_FILES)
            self.dnd_bind("<<Drop>>", self._on_drop)
            self.dnd_bind("<<DropEnter>>", lambda _e: self._tint(HOT))
            self.dnd_bind("<<DropLeave>>", lambda _e: self._tint(PANEL))

    def _tint(self, colour):
        for w in (self, self.lbl_title, self.lbl_path, self.lbl_info):
            w.configure(bg=colour)

    def _on_drop(self, event):
        self._tint(PANEL)
        paths = [p for p in parse_drop(event.data) if os.path.isfile(p)]
        if paths:
            self.on_files(self, paths)

    def browse(self):
        paths = filedialog.askopenfilenames(
            title=f"Choose the {self.title} file",
            filetypes=[("xyz / text", "*.xyz *.txt *.dat *.d *"), ("all", "*")])
        if paths:
            self.on_files(self, list(paths))

    def show(self, path, info, colour="#374151"):
        self.lbl_path.configure(text=os.path.basename(path), fg="#111827")
        self.lbl_info.configure(text=info, fg=colour)


# ============================================================================
# Main window
# ============================================================================
class App:
    def __init__(self, root):
        self.root = root
        root.title("Tomographic map merger")
        root.configure(bg=BG)
        root.minsize(900, 700)

        self.data = {"west": None, "east": None}     # (lon, lat, vel)
        self.path = {"west": "", "east": ""}
        self.q = queue.Queue()
        self._sync = False
        self.last_qc = None
        self.busy = False

        self.varL = tk.DoubleVar(value=1.0)
        self.varLtxt = tk.StringVar(value="1.00")
        self.varCW = tk.StringVar(value="")
        self.varCE = tk.StringVar(value="")
        self.varOut = tk.StringVar(value="")
        self.varQC = tk.BooleanVar(value=True)
        self.varDiag = tk.BooleanVar(value=False)

        self._build()
        self.varL.trace_add("write", self._on_L_var)
        for v in (self.varCW, self.varCE):
            v.trace_add("write", lambda *_a: self._refresh())
        self.root.after(80, self._pump)

    # ---------------------------------------------------------------- layout
    def _build(self):
        r = self.root

        tk.Label(r, text="Merge two overlapping tomographic maps", bg=BG,
                 fg=ACCENT, font=("Segoe UI", 15, "bold")).pack(
                     anchor="w", padx=16, pady=(14, 0))
        tk.Label(r, text="Gaussian distance-weighted blend in the overlap zone; "
                         "nodes seen by only one array pass through untouched.",
                 bg=BG, fg="#4b5563", font=("Segoe UI", 9)).pack(
                     anchor="w", padx=16, pady=(1, 10))

        # --- inputs
        row = tk.Frame(r, bg=BG); row.pack(fill="x", padx=16)
        self.panelW = DropPanel(row, "WEST  .xyz", self._on_files)
        self.panelE = DropPanel(row, "EAST  .xyz", self._on_files)
        self.panelW.pack(side="left", fill="both", expand=True, padx=(0, 6))
        self.panelE.pack(side="left", fill="both", expand=True, padx=(6, 0))

        self.lbl_overlap = tk.Label(r, text="waiting for both files", bg=BG,
                                    fg="#6b7280", font=("Segoe UI", 9))
        self.lbl_overlap.pack(anchor="w", padx=16, pady=(8, 0))

        # --- parameters
        box = tk.LabelFrame(r, text=" blending ", bg=BG, fg=ACCENT,
                            font=("Segoe UI", 10, "bold"), bd=1,
                            relief="groove")
        box.pack(fill="x", padx=16, pady=10)

        cr = tk.Frame(box, bg=BG); cr.pack(fill="x", padx=10, pady=(8, 2))
        tk.Label(cr, text="array centre lon    west", bg=BG,
                 font=("Segoe UI", 9)).pack(side="left")
        ttk.Entry(cr, textvariable=self.varCW, width=9).pack(side="left", padx=6)
        tk.Label(cr, text="east", bg=BG, font=("Segoe UI", 9)).pack(side="left")
        ttk.Entry(cr, textvariable=self.varCE, width=9).pack(side="left", padx=6)
        tk.Label(cr, text="(defaults: mid-longitude of each map)", bg=BG,
                 fg="#6b7280", font=("Segoe UI", 8)).pack(side="left", padx=6)

        sr = tk.Frame(box, bg=BG); sr.pack(fill="x", padx=10, pady=(6, 2))
        tk.Label(sr, text="smoothing length  L", bg=BG,
                 font=("Segoe UI", 9)).pack(side="left")
        self.scale = ttk.Scale(sr, from_=0.05, to=10.0, orient="horizontal",
                               variable=self.varL)
        self.scale.pack(side="left", fill="x", expand=True, padx=8)
        e = ttk.Entry(sr, textvariable=self.varLtxt, width=8)
        e.pack(side="left")
        e.bind("<Return>", self._on_L_text)
        e.bind("<FocusOut>", self._on_L_text)
        tk.Label(sr, text="deg", bg=BG, font=("Segoe UI", 9)).pack(side="left",
                                                                   padx=(3, 8))
        ttk.Button(sr, text="Recommend", command=self._recommend).pack(side="left")

        self.lbl_kernel = tk.Label(box, text="", bg=BG, fg="#374151",
                                   font=("Consolas", 9), justify="left")
        self.lbl_kernel.pack(anchor="w", padx=10, pady=(4, 2))
        self.lbl_verdict = tk.Label(box, text="", bg=BG, fg="#374151",
                                    font=("Segoe UI", 9), justify="left",
                                    wraplength=820)
        self.lbl_verdict.pack(anchor="w", padx=10, pady=(0, 8))

        # --- live preview
        if HAVE_MPL:
            self.fig = Figure(figsize=(7.4, 2.0), dpi=100)
            self.ax = self.fig.add_subplot(111)
            self.fig.subplots_adjust(left=.07, right=.99, top=.88, bottom=.26)
            self.canvas = FigureCanvasTkAgg(self.fig, master=r)
            self.canvas.get_tk_widget().pack(fill="x", padx=16)
        else:
            self.fig = self.ax = self.canvas = None

        # --- output
        orow = tk.Frame(r, bg=BG); orow.pack(fill="x", padx=16, pady=(10, 0))
        tk.Label(orow, text="output", bg=BG, font=("Segoe UI", 9)).pack(side="left")
        ttk.Entry(orow, textvariable=self.varOut).pack(side="left", fill="x",
                                                       expand=True, padx=6)
        ttk.Button(orow, text="Save as...", command=self._save_as).pack(side="left")

        crow = tk.Frame(r, bg=BG); crow.pack(fill="x", padx=16, pady=(4, 0))
        ttk.Checkbutton(crow, text="QC figure (.png)",
                        variable=self.varQC).pack(side="left")
        ttk.Checkbutton(crow, text="per-node weights (_diag.txt)",
                        variable=self.varDiag).pack(side="left", padx=12)

        brow = tk.Frame(r, bg=BG); brow.pack(fill="x", padx=16, pady=(10, 0))
        self.btn = tk.Button(brow, text="MERGE", command=self._merge,
                             bg=ACCENT, fg="white", font=("Segoe UI", 11, "bold"),
                             relief="flat", padx=26, pady=7,
                             activebackground="#16386a", activeforeground="white",
                             disabledforeground="#c7d2e4")
        self.btn.pack(side="left")
        self.btn.configure(state="disabled")
        self.btn_qc = ttk.Button(brow, text="Open QC figure",
                                 command=self._open_qc, state="disabled")
        self.btn_qc.pack(side="left", padx=10)
        self.prog = ttk.Progressbar(brow, mode="indeterminate", length=170)
        self.prog.pack(side="left", padx=10)

        lf = tk.Frame(r, bg=BG); lf.pack(fill="both", expand=True, padx=16,
                                         pady=(10, 14))
        self.log = tk.Text(lf, height=9, bg="#0f172a", fg="#e2e8f0",
                           font=("Consolas", 9), relief="flat", wrap="none",
                           insertbackground="#e2e8f0")
        sb = ttk.Scrollbar(lf, orient="vertical", command=self.log.yview)
        self.log.configure(yscrollcommand=sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        sb.pack(side="left", fill="y")
        self.log.tag_configure("err", foreground="#fca5a5")
        self.log.tag_configure("ok", foreground="#86efac")

        if not HAVE_DND:
            self._say("drag-and-drop is off (pip install tkinterdnd2 to enable) - "
                      "click a panel to browse instead")
        self._refresh()

    # ------------------------------------------------------------- log / bg
    def _say(self, text, tag=None):
        self.log.insert("end", text.rstrip() + "\n", tag or ())
        self.log.see("end")

    def _pump(self):
        """Drain messages posted by worker threads onto the Tk thread."""
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "log":
                    self._say(payload)
                elif kind == "loaded":
                    self._apply_loaded(*payload)
                elif kind == "merged":
                    self._merge_done(payload)
                elif kind == "error":
                    self._stop_busy()
                    self._say(payload, "err")
                    messagebox.showerror("Error", payload)
        except queue.Empty:
            pass
        self.root.after(80, self._pump)

    def _bg(self, fn):
        threading.Thread(target=fn, daemon=True).start()

    def _start_busy(self):
        self.busy = True
        self.btn.configure(state="disabled", text="working...")
        self.prog.start(12)

    def _stop_busy(self):
        self.busy = False
        self.prog.stop()
        self.btn.configure(text="MERGE")
        self._refresh()

    # --------------------------------------------------------------- inputs
    def _on_files(self, panel, paths):
        side = "west" if panel is self.panelW else "east"

        def work():
            try:
                loaded = [(p, mt.read_xyz(p)) for p in paths[:2]]
            except Exception as exc:
                self.q.put(("error", str(exc))); return
            if len(loaded) == 2:
                # two files at once: sort them into the right panels
                w, e = order_by_longitude(loaded)
                self.q.put(("loaded", ("west", w[0], w[1])))
                self.q.put(("loaded", ("east", e[0], e[1])))
            else:
                self.q.put(("loaded", (side, loaded[0][0], loaded[0][1])))

        self._bg(work)

    def _apply_loaded(self, side, path, data):
        self.data[side] = data
        self.path[side] = path
        panel = self.panelW if side == "west" else self.panelE
        panel.show(path, summarise(*data))
        self._say(f"loaded {side}: {path}  ({len(data[0])} nodes)")

        if not self.varOut.get():
            self.varOut.set(default_output(self.path["west"] or path))

        if self.data["west"] is not None and self.data["east"] is not None:
            cW, cE = default_centres(self.data["west"], self.data["east"])
            self._sync = True
            self.varCW.set(f"{cW}"); self.varCE.set(f"{cE}")
            self._sync = False
            self._recommend()
        self._refresh()

    def _centres(self):
        try:
            return float(self.varCW.get()), float(self.varCE.get())
        except ValueError:
            return None, None

    def _ready(self):
        return (self.data["west"] is not None and self.data["east"] is not None
                and self._centres()[0] is not None)

    # ------------------------------------------------------------ L control
    def _on_L_var(self, *_a):
        if self._sync:
            return
        self._sync = True
        self.varLtxt.set(f"{self.varL.get():.2f}")
        self._sync = False
        self._refresh()

    def _on_L_text(self, *_a):
        try:
            v = float(self.varLtxt.get())
        except ValueError:
            self.varLtxt.set(f"{self.varL.get():.2f}"); return
        if v <= 0:
            messagebox.showwarning("Smoothing length",
                                   "L must be greater than zero.")
            self.varLtxt.set(f"{self.varL.get():.2f}"); return
        lo, hi = float(self.scale.cget("from")), float(self.scale.cget("to"))
        if v > hi:                       # let a typed value widen the slider
            self.scale.configure(to=round(v * 1.2, 2))
        if v < lo:
            self.scale.configure(from_=max(0.01, v * 0.8))
        self.varL.set(v)

    def _recommend(self):
        if not self._ready():
            return
        cW, cE = self._centres()
        ov_lo, ov_hi = overlap_span(self.data["west"], self.data["east"])
        lo, hi = slider_range(cW, cE, ov_lo, ov_hi)
        self.scale.configure(from_=lo, to=hi)
        Ls = mt.suggest_L(cW, cE, ov_lo, ov_hi)
        if Ls is None:
            self._say("! the midpoint between the centres is outside the overlap; "
                      "no L keeps both edges clean - check the centres", "err")
            return
        self.varL.set(float(np.floor(Ls * 100) / 100))   # floor, never round up

    # ------------------------------------------------------------- refresh
    def _refresh(self):
        ready = self._ready()
        if not self.busy:
            self.btn.configure(state="normal" if ready else "disabled")
        if not ready:
            self.lbl_kernel.configure(text="")
            self.lbl_verdict.configure(text="")
            if self.ax:
                self.ax.clear(); self.canvas.draw_idle()
            return

        cW, cE = self._centres()
        ov_lo, ov_hi = overlap_span(self.data["west"], self.data["east"])
        L = max(self.varL.get(), 1e-6)
        s = mt.handover_scale(L, cW, cE)

        if ov_hi <= ov_lo:
            self.lbl_overlap.configure(
                text="the two maps do not overlap in longitude - "
                     "they will simply be concatenated", fg=WARN_C)
        else:
            self.lbl_overlap.configure(
                text=f"longitude overlap  {ov_lo:.2f} .. {ov_hi:.2f}   "
                     f"({ov_hi-ov_lo:.2f} deg wide)", fg="#374151")

        wlo, wehi, verdict = edge_quality(L, cW, cE, ov_lo, ov_hi)
        self.lbl_kernel.configure(
            text=f"hand-over scale  s = L^2/(2*{abs(cE-cW):.2f}) = {s:.2f} deg     "
                 f"10-90% width {4.394*s:.2f} deg")

        if verdict == "none":
            self.lbl_verdict.configure(text="", fg="#374151")
        elif verdict == "ok":
            self.lbl_verdict.configure(
                text=f"edges clean:  {wlo:.0%} west at lon {ov_lo:.2f},  "
                     f"{wehi:.0%} east at lon {ov_hi:.2f}", fg=OK_C)
        else:
            self.lbl_verdict.configure(
                text=f"too smooth for this overlap: only {wlo:.0%} west at lon "
                     f"{ov_lo:.2f} and {wehi:.0%} east at lon {ov_hi:.2f}, so the "
                     f"merged map will step where the overlap ends. "
                     f"Press Recommend for the largest safe L.", fg=BAD_C)

        self._draw_preview(L, cW, cE, ov_lo, ov_hi)

    def _draw_preview(self, L, cW, cE, ov_lo, ov_hi):
        if not self.ax:
            return
        lonmin = min(self.data["west"][0].min(), self.data["east"][0].min())
        lonmax = max(self.data["west"][0].max(), self.data["east"][0].max())
        xs = np.linspace(lonmin, lonmax, 800)
        w = mt.west_weight(xs, L, cW, cE)
        ax = self.ax
        ax.clear()
        if ov_hi > ov_lo:
            ax.axvspan(ov_lo, ov_hi, color="0.5", alpha=.18, label="overlap")
        ax.plot(xs, w, lw=2, label="west")
        ax.plot(xs, 1 - w, lw=2, label="east")
        ax.axhline(mt.EDGE_TOL, ls="--", lw=.8, c="0.5")
        ax.axhline(1 - mt.EDGE_TOL, ls="--", lw=.8, c="0.5")
        ax.axvline(cW, ls=":", lw=1, c="C0"); ax.axvline(cE, ls=":", lw=1, c="C1")
        ax.set_xlim(lonmin, lonmax); ax.set_ylim(-.05, 1.05)
        ax.set_xlabel("longitude", fontsize=8)
        ax.set_ylabel("weight", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="center left", ncol=3, framealpha=.9)
        ax.grid(alpha=.25)
        self.canvas.draw_idle()

    # --------------------------------------------------------------- output
    def _save_as(self):
        p = filedialog.asksaveasfilename(defaultextension=".xyz",
                                         initialfile="merged.xyz",
                                         filetypes=[("xyz", "*.xyz"),
                                                    ("all", "*")])
        if p:
            self.varOut.set(p)

    def _open_qc(self):
        if self.last_qc and os.path.isfile(self.last_qc):
            if not open_externally(self.last_qc):
                messagebox.showinfo("QC figure", self.last_qc)

    # ---------------------------------------------------------------- merge
    def _merge(self):
        if not self._ready() or self.busy:
            return
        out = self.varOut.get().strip()
        if not out:
            messagebox.showwarning("Output", "Choose an output file first.")
            return
        if os.path.isfile(out) and not messagebox.askyesno(
                "Overwrite?", f"{os.path.basename(out)} already exists.\n"
                              f"Overwrite it?"):
            return

        cW, cE = self._centres()
        L = self.varL.get()
        west, east = self.data["west"], self.data["east"]
        want_qc, want_diag = self.varQC.get(), self.varDiag.get()

        self._start_busy()
        self._say(f"\nmerging with L = {L:.2f} deg, centres {cW} / {cE} ...")

        def work():
            try:
                import io, contextlib
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    lon, lat, vel, wW, src = mt.gaussian_merge(west, east, L,
                                                               cW, cE)
                    ov, diff = mt.overlap_diff(west, east, lon, lat, src)
                    with open(out, "w") as fh:
                        for x, y, v in zip(lon, lat, vel):
                            fh.write(f"{x:12.6f} {y:12.7f} {v:14.8f}\n")
                    diag = None
                    if want_diag:
                        diag = os.path.splitext(out)[0] + "_diag.txt"
                        with open(diag, "w") as fh:
                            fh.write("# lon lat vel w_W source(W/E/B)\n")
                            for x, y, v, w_, s_ in zip(lon, lat, vel, wW, src):
                                fh.write(f"{x:12.6f} {y:12.7f} {v:14.8f} "
                                         f"{w_:9.6f} {s_}\n")
                    png = None
                    if want_qc:
                        png = os.path.splitext(out)[0] + "_qc.png"
                        mt.qc_plot(west, east, (lon, lat, vel), wW, src, diff,
                                   ov, L, cW, cE, png)
                    sr = mt.seam_report(lon, wW, ov, diff) if ov.any() else None
                notes = buf.getvalue().strip()
                self.q.put(("merged", dict(
                    out=out, diag=diag, png=png, n=len(lon),
                    nW=int((src == "W").sum()), nE=int((src == "E").sum()),
                    nB=int(ov.sum()), diff=diff, seam=sr, notes=notes)))
            except Exception:
                self.q.put(("error", traceback.format_exc(limit=3)))

        self._bg(work)

    def _merge_done(self, r):
        self._stop_busy()
        if r["notes"]:
            self._say(r["notes"])
        self._say(f"{r['n']} nodes written:  {r['nW']} west-only, "
                  f"{r['nE']} east-only, {r['nB']} blended")
        if r["nB"]:
            d = r["diff"]
            self._say(f"overlap disagreement (east - west):  mean {d.mean():+.4f}"
                      f"   rms {np.sqrt((d**2).mean()):.4f}"
                      f"   range {d.min():+.4f} .. {d.max():+.4f}")
            s = r["seam"]
            if s:
                self._say(f"residual step at overlap edges:  "
                          f"lon {s['lo_edge']:.2f} rms {s['rms_w']:.4f} "
                          f"(max {s['max_w']:.4f})   "
                          f"lon {s['hi_edge']:.2f} rms {s['rms_e']:.4f} "
                          f"(max {s['max_e']:.4f})")
        self._say(f"wrote {r['out']}", "ok")
        if r["diag"]:
            self._say(f"wrote {r['diag']}", "ok")
        if r["png"]:
            self._say(f"wrote {r['png']}", "ok")
            self.last_qc = r["png"]
            self.btn_qc.configure(state="normal")


def main():
    root = TkinterDnD.Tk() if HAVE_DND else tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
