#!/usr/bin/env python3
"""Print the LaTeX text of boxed statements by their index number.

Usage: python3 -I boxes.py <compendium_dir> <number> [<number> ...]
Numbering follows the source's index (sections counted continuously, lettered after
\\appendix, I.n among the instruments; the five boxed environments share one counter
per section). Read-only; nothing in the source is executed.
"""
import re, sys, pathlib

ENV = "chtheorem|chconditional|chwall|chdefinition|chtable"

def balanced(s, i, o="{", c="}"):
    d = 0
    for j in range(i, len(s)):
        if s[j] == o: d += 1
        elif s[j] == c:
            d -= 1
            if d == 0: return s[i+1:j], j+1
    return s[i+1:], len(s)

def boxes(texdir):
    order = ["front"] + [p.stem for p in sorted((texdir/"parts").glob("*.tex"))] + ["back", "instruments"]
    tok = re.compile(r"\\part\*\{|\\appendix\b|\\section\{|\\begin\{(" + ENV + r")\}")
    secnum, appendix, instr, cur, k = 0, False, False, None, 0
    out = {}
    for f in order:
        p = texdir/(f + ".tex") if f in ("front", "back", "instruments") else texdir/"parts"/(f + ".tex")
        s = p.read_text(encoding="utf-8"); pos = 0
        while True:
            m = tok.search(s, pos)
            if not m: break
            t = m.group(0)
            if t.startswith("\\part*"):
                title, e = balanced(s, m.end()-1)
                if title.strip() == "Instruments": instr, appendix, secnum = True, False, 0
                pos = e
            elif t.startswith("\\appendix"):
                appendix, secnum, pos = True, 0, m.end()
            elif t.startswith("\\section{"):
                title, e = balanced(s, m.end()-1); secnum += 1; k = 0
                cur = f"I.{secnum}" if instr else (chr(64+secnum) if appendix else str(secnum))
                pos = e
            else:
                env = m.group(1); k += 1; e = m.end(); opt = ""
                if s[e:e+1] == "[": opt, e = balanced(s, e, "[", "]")
                end = s.find("\\end{" + env + "}", e)
                out[f"{cur}.{k}"] = dict(env=env, title=opt, body=s[e:end].strip(), file=p.name,
                                         line=s.count("\n", 0, m.start()) + 1)
                pos = end + 1
    return out

if __name__ == "__main__":
    B = boxes(pathlib.Path(sys.argv[1])/"tex")
    for n in sys.argv[2:]:
        b = B.get(n)
        if not b: print(f"== {n}: not found"); continue
        print(f"== {n} [{b['env']}] {b['title']}  ({b['file']}:{b['line']})\n{b['body']}\n")
