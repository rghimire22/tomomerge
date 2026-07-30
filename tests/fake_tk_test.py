"""
Headless test rig for merge_tomo_gui.py.

There is no tkinter in this sandbox, so a stand-in is injected into
sys.modules.  It mimics the parts of the Tk API the GUI actually touches -
including Tk's rule that a trailing underscore is stripped from option names,
so cget("from") sees what configure(from_=...) wrote.  That lets every
callback, the worker function and the queue protocol run for real; only Tk's
own rendering is out of reach.
"""
import sys, types, os

# --------------------------------------------------------------- fake tkinter
class _W:
    def __init__(self, master=None, **kw):
        self.master = master
        self._kw = {}
        self.children = []
        self.configure(**kw)
        if isinstance(master, _W):
            master.children.append(self)

    def configure(self, **kw):
        for k, v in kw.items():
            self._kw[k[:-1] if k.endswith("_") else k] = v
    config = configure

    def cget(self, k):
        return self._kw[k]

    def __getitem__(self, k):
        return self._kw[k]

    def pack(self, **k): pass
    def grid(self, **k): pass
    def place(self, **k): pass
    def pack_forget(self): pass
    def bind(self, *a, **k): pass
    def focus_set(self): pass
    def update_idletasks(self): pass
    def winfo_exists(self): return 1
    def set(self, *a): pass          # Scrollbar.set
    def get(self, *a): return 0.0
    # tkinterdnd2 surface (never reached: HAVE_DND is False here)
    def drop_target_register(self, *a): pass
    def dnd_bind(self, *a): pass


class _Text(_W):
    def __init__(self, *a, **k):
        super().__init__(*a, **k); self.buf = []
    def insert(self, where, text, tags=()): self.buf.append(text)
    def see(self, *a): pass
    def tag_configure(self, *a, **k): pass
    def yview(self, *a): pass
    def text(self): return "".join(self.buf)


class _Prog(_W):
    def start(self, *a): self._kw["_running"] = True
    def stop(self, *a): self._kw["_running"] = False


class _Var:
    def __init__(self, master=None, value=None):
        self._v = value; self._cbs = []
    def get(self): return self._v
    def set(self, v):
        self._v = v
        for cb in list(self._cbs):
            cb("", "", "write")
    def trace_add(self, mode, cb): self._cbs.append(cb)


class _Root(_W):
    def __init__(self, *a, **k):
        super().__init__(None, **k); self.pending = []
    def title(self, *a): pass
    def minsize(self, *a): pass
    def mainloop(self): pass
    def after(self, ms, fn=None, *a):
        self.pending.append(fn)          # never auto-fires; test drives _pump
        return "id"


tk = types.ModuleType("tkinter")
for name in ("Frame", "Label", "Button", "LabelFrame", "Canvas", "Widget",
             "BaseWidget", "Toplevel"):
    setattr(tk, name, type(name, (_W,), {}))
tk.Tk = _Root
tk.Text = _Text
tk.StringVar = lambda master=None, value="": _Var(master, value)
tk.DoubleVar = lambda master=None, value=0.0: _Var(master, value)
tk.BooleanVar = lambda master=None, value=False: _Var(master, value)

ttk = types.ModuleType("tkinter.ttk")
for name in ("Entry", "Button", "Scale", "Checkbutton", "Scrollbar",
             "Label", "Frame", "Combobox"):
    setattr(ttk, name, type(name, (_W,), {}))
ttk.Progressbar = _Prog

CALLS = {"info": [], "warn": [], "error": [], "askyesno": []}
messagebox = types.ModuleType("tkinter.messagebox")
messagebox.showinfo = lambda t, m, **k: CALLS["info"].append((t, m))
messagebox.showwarning = lambda t, m, **k: CALLS["warn"].append((t, m))
messagebox.showerror = lambda t, m, **k: CALLS["error"].append((t, m))
messagebox.askyesno = lambda t, m, **k: (CALLS["askyesno"].append((t, m)), True)[1]

filedialog = types.ModuleType("tkinter.filedialog")
filedialog.askopenfilenames = lambda **k: ()
filedialog.asksaveasfilename = lambda **k: ""

tk.ttk, tk.messagebox, tk.filedialog = ttk, messagebox, filedialog
sys.modules["tkinter"] = tk
sys.modules["tkinter.ttk"] = ttk
sys.modules["tkinter.messagebox"] = messagebox
sys.modules["tkinter.filedialog"] = filedialog

sys.path.insert(0, "/sessions/trusting-charming-babbage/mnt/Merging")
import numpy as np
import merge_tomo as mt
import merge_tomo_gui as g

OUT = "/sessions/trusting-charming-babbage/mnt/outputs"
W = os.path.join(OUT, "west.xyz")
E = os.path.join(OUT, "east.xyz")
fails = []


def check(name, cond, extra=""):
    print(("  PASS  " if cond else "  FAIL  ") + name + (f"   {extra}" if extra else ""))
    if not cond:
        fails.append(name)


print("\n--- pure helpers ---")
check("parse_drop: plain", g.parse_drop("C:/a.xyz C:/b.xyz") ==
      ["C:/a.xyz", "C:/b.xyz"])
check("parse_drop: braces + spaces",
      g.parse_drop("{C:/My Data/west.xyz} C:/east.xyz") ==
      ["C:/My Data/west.xyz", "C:/east.xyz"])
check("parse_drop: single braced",
      g.parse_drop("{/home/rg/vel 02.xyz}") == ["/home/rg/vel 02.xyz"])
check("parse_drop: empty", g.parse_drop("") == [])

w, e = mt.read_xyz(W), mt.read_xyz(E)
check("order_by_longitude", [p for p, _ in g.order_by_longitude(
    [(E, e), (W, w)])] == [W, E])
cW, cE = g.default_centres(w, e)
check("default_centres", (cW, cE) == (-119.0, -110.0), f"{cW} {cE}")
ov = g.overlap_span(w, e)
check("overlap_span", ov == (-117.0, -112.0), str(ov))

Ls = mt.suggest_L(cW, cE, *ov)
wl, we_, verdict = g.edge_quality(Ls, cW, cE, *ov)
check("edge_quality at suggested L is 'ok'", verdict == "ok",
      f"w_W={wl:.4f} w_E={we_:.4f}")
check("edge_quality flags an over-smooth L",
      g.edge_quality(Ls * 2, cW, cE, *ov)[2] == "loose")
check("edge_quality: no overlap -> 'none'",
      g.edge_quality(4.0, cW, cE, -100.0, -105.0)[2] == "none")
lo, hi = g.slider_range(cW, cE, *ov)
check("slider_range brackets the recommendation", lo < Ls < hi, f"{lo}..{hi}")

print("\n--- app wiring ---")
root = tk.Tk()
app = g.App(root)
app._bg = lambda fn: fn()          # run 'threads' inline for determinism

app._on_files(app.panelW, [E, W])  # two files at once, deliberately wrong order
app._pump()
check("two dropped files are sorted into the right panels",
      app.path["west"] == W and app.path["east"] == E,
      f"west={os.path.basename(app.path['west'])}")
check("centres auto-filled", (app.varCW.get(), app.varCE.get()) ==
      ("-119.0", "-110.0"))
check("L auto-set to floored recommendation",
      abs(app.varL.get() - np.floor(Ls * 100) / 100) < 1e-12,
      f"L={app.varL.get()}")
check("entry text tracks the slider", app.varLtxt.get() ==
      f"{app.varL.get():.2f}")
check("verdict is green at the default L",
      app.lbl_verdict.cget("fg") == g.OK_C)
check("MERGE enabled once both files are in",
      app.btn.cget("state") == "normal")

app.varL.set(Ls * 3)
check("verdict turns red when over-smoothed",
      app.lbl_verdict.cget("fg") == g.BAD_C)

app.varLtxt.set("-2"); app._on_L_text()
check("negative L rejected and reverted", len(CALLS["warn"]) == 1 and
      abs(app.varL.get() - Ls * 3) < 1e-9)
app.varLtxt.set("40"); app._on_L_text()
check("typing past the slider max widens the slider",
      float(app.scale.cget("to")) >= 40 and abs(app.varL.get() - 40) < 1e-9)

app.varL.set(float(np.floor(Ls * 100) / 100))

print("\n--- merge through the GUI equals the CLI ---")
gui_out = os.path.join(OUT, "gui_merged.xyz")
app.varOut.set(gui_out)
app.varQC._v = False               # skip the slow figure in this check
app.varDiag._v = True
app._merge()
app._pump()
check("no error dialog", not CALLS["error"], str(CALLS["error"])[:200])
check("output written", os.path.isfile(gui_out))
check("diag written", os.path.isfile(os.path.join(OUT, "gui_merged_diag.txt")))

cli = np.loadtxt(os.path.join(OUT, "merged.xyz"))
gui = np.loadtxt(gui_out)
check("GUI output is byte-identical to the CLI output",
      cli.shape == gui.shape and np.array_equal(cli, gui))

log = app.log.text()
check("log reports node counts", "34881 nodes written" in log)
check("log reports the seam", "residual step at overlap edges" in log)
check("Open-QC button stays disabled when QC is off",
      app.btn_qc.cget("state") == "disabled")

print("\n--- non-overlapping maps must not crash ---")
np.savetxt(os.path.join(OUT, "far.xyz"),
           np.c_[np.repeat(np.arange(-90, -85, .5), 5),
                 np.tile(np.arange(40, 42.5, .5), 10),
                 np.full(50, 4.2)])
app2 = g.App(tk.Tk()); app2._bg = lambda fn: fn()
app2._on_files(app2.panelW, [W]); app2._pump()
app2._on_files(app2.panelE, [os.path.join(OUT, "far.xyz")]); app2._pump()
check("non-overlap is reported, app still usable",
      "do not overlap" in app2.lbl_overlap.cget("text"))
app2.varOut.set(os.path.join(OUT, "far_merged.xyz"))
app2.varQC._v = False
app2._merge(); app2._pump()
check("concatenation still succeeds with zero blended nodes",
      os.path.isfile(os.path.join(OUT, "far_merged.xyz")) and
      "0 blended" in app2.log.text())

print("\n" + ("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}"))
sys.exit(1 if fails else 0)
