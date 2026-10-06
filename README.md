# Sky Sketch

Connect stars in a night sky. When you finish, it shows the closest of the 88
modern constellations and how other cultures drew the same stars.

Part of my design playground — see it live at
[portfolio.lakshmisrinath.com](https://portfolio.lakshmisrinath.com/) (Playground → Constellation).

Constellation figures and cross-cultural names come from the
[Stellarium sky cultures](https://github.com/Stellarium/stellarium-skycultures) project.

## Run it

Static page, no build step, but it loads data with `fetch`, so serve it over HTTP:

```bash
python -m http.server 8000
# then open http://localhost:8000
```
