# frontend

A minimal static viewer for the variant × therapy table. Next.js with `output: export`, so it
builds to plain files and deploys to Cloudflare Pages without a server. It reads one file,
`public/riborescue.json`, written by `riborescue export-web` — the app never touches a BAM or a
results table directly, only that compact JSON.

The table (sort, filter, search) is TanStack Table; the app owns only the rendering and the two
status banners, so there is little here to maintain and it can be replaced wholesale.

```
pixi run export-web-example   # regenerate public/riborescue.json from results/
pixi run app-install          # npm install, once
pixi run app-dev              # dev server at http://localhost:3000
pixi run app-build           # static export to frontend/out/
```

To point it at a real dataset rather than the committed example, run `riborescue export-web` without
`--sample` over the full `results/landscape.tsv` and drop the JSON in `public/`.
