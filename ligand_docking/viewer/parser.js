/* Minimal AutoDock PDBQT helpers for the docking viewer.
   Works in the browser (window.DockParser) and in node (module.exports). */
(function (root) {
  "use strict";

  const ELEMENT_BY_TYPE = {
    A: "C", C: "C", N: "N", NA: "N", NS: "N", OA: "O", O: "O",
    SA: "S", S: "S", HD: "H", H: "H", HS: "H", P: "P", F: "F",
    CL: "Cl", BR: "Br", I: "I", SI: "Si", MG: "Mg", ZN: "Zn",
    FE: "Fe", CA: "Ca", MN: "Mn", SE: "Se", NA_: "Na",
  };
  const KNOWN = ["Cl", "Br", "Si", "Na", "Mg", "Zn", "Fe", "Mn", "Ca", "Se",
                 "F", "I", "P", "S", "O", "N", "C", "H"];

  function elementOf(atomType, atomName) {
    const key = (atomType || "").trim().toUpperCase();
    if (ELEMENT_BY_TYPE[key]) return ELEMENT_BY_TYPE[key];
    for (const symbol of KNOWN) {
      if (key.startsWith(symbol.toUpperCase())) return symbol;
    }
    const first = (atomName || "").trim().charAt(0).toUpperCase();
    return first || "C";
  }

  function atomLineToPdb(line) {
    const record = line.slice(0, 6);
    const name = line.slice(12, 16);
    const type = line.length >= 79 ? line.slice(77, 79) : line.slice(77);
    const element = elementOf(type, name);
    const body = line.slice(0, 54).padEnd(54, " ");
    const occupancy = (line.slice(54, 60) || "").trim() ? line.slice(54, 60) : "  1.00";
    const bfactor = (line.slice(60, 66) || "").trim() ? line.slice(60, 66) : "  0.00";
    const head = record + body.slice(6, 54) + occupancy + bfactor;   // 66 columns
    const label = element.length === 1
      ? element.toUpperCase()
      : element.charAt(0).toUpperCase() + element.slice(1).toLowerCase();
    return head + " ".repeat(10) + label.padStart(2, " ");
  }

  /** Convert any PDBQT text (receptor or ligand, single or multi model) to PDB. */
  function pdbqtToPdb(text) {
    const out = [];
    for (const raw of text.split(/\r?\n/)) {
      if (raw.startsWith("ATOM") || raw.startsWith("HETATM")) out.push(atomLineToPdb(raw));
    }
    return out.join("\n") + "\nEND\n";
  }

  /** Parse a Vina `--out` file into poses with scores. */
  function parseVinaOut(text) {
    const poses = [];
    let current = null;
    let remark = null;
    const flush = () => {
      if (current && current.length) {
        poses.push({
          mode: poses.length + 1,
          affinity: remark ? remark.affinity : null,
          rmsd_lb: remark ? remark.rmsd_lb : null,
          rmsd_ub: remark ? remark.rmsd_ub : null,
          pdb: current.join("\n") + "\nEND\n",
        });
      }
      current = null;
      remark = null;
    };

    for (const raw of text.split(/\r?\n/)) {
      if (raw.startsWith("MODEL")) {
        current = [];
        remark = null;
      } else if (raw.startsWith("REMARK VINA RESULT:")) {
        const f = raw.trim().split(/\s+/);
        remark = { affinity: parseFloat(f[3]), rmsd_lb: parseFloat(f[4]), rmsd_ub: parseFloat(f[5]) };
      } else if (raw.startsWith("ATOM") || raw.startsWith("HETATM")) {
        if (!current) current = [];
        current.push(atomLineToPdb(raw));
      } else if (raw.startsWith("ENDMDL")) {
        flush();
      }
    }
    flush();
    return poses;
  }

  /** Guess what a dropped file is from its name and first lines. */
  function sniff(name, text) {
    const lower = (name || "").toLowerCase();
    const head = text.slice(0, 4000);
    if (lower.endsWith(".pdbqt") || /^TORSDOF/m.test(head) || /^ROOT\b/m.test(head)) {
      return /^REMARK VINA RESULT/m.test(text) ? "poses" : "pdbqt";
    }
    if (lower.endsWith(".sdf") || lower.endsWith(".mol") || text.includes("$$$$")) return "sdf";
    if (lower.endsWith(".mol2")) return "mol2";
    if (lower.endsWith(".cif") || lower.endsWith(".mmcif")) return "cif";
    if (lower.endsWith(".pdb") || lower.endsWith(".ent") || /^(ATOM|HETATM)/m.test(text)) return "pdb";
    return "unknown";
  }

  const api = { elementOf, atomLineToPdb, pdbqtToPdb, parseVinaOut, sniff };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.DockParser = api;
})(typeof window !== "undefined" ? window : globalThis);
