from pathlib import Path

p = Path("dmfb_checker.py")

txt = p.read_text()

old = '''da = new_id("S")
                db = new_id("S")
                new_positions[da] = (r1, c1)
                new_positions[db] = (r2, c2)'''

new = '''dm = new_id("M")
                new_positions[dm] = (r1, c1)'''

txt = txt.replace(old, new)

Path("dmfb_checker_v2.py").write_text(txt)

print("✅ Wrote dmfb_checker_v2.py")