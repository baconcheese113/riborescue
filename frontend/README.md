# frontend

A minimal static viewer for the variant × therapy table. Next.js with `output: export`, so it
builds to plain files and deploys to Cloudflare Pages without a server. It reads one file,
`public/riborescue.json`, written by `riborescue export-web` — the app never touches a BAM or a
results table directly, only that compact JSON.

Pages builds it from this directory: root `frontend`, build command `npm run build`, output
directory `out`. Both viewer JSONs are committed, so the build needs neither Pixi nor the pipeline.

The `overrides` in `package.json` are not optional. Next pins `postcss` at exactly 8.4.31 and
`sharp` at `^0.34.5`, and both carry advisories in every published Next release, so the versions
have to be forced from here.

The app owns the visual explanations and lookup interactions. Scientific values stay in the
exported JSON rather than the components.

```
pixi run export-web-example   # regenerate public/riborescue.json from results/
pixi run app-install          # npm install, once
pixi run app-dev              # dev server at http://localhost:3000
pixi run app-build            # static export to frontend/out/
pixi run app-preview          # build, then serve the export at http://localhost:3001
```

`app-preview` serves the same static export that deploys. Use it for browser validation.

To point it at a real dataset rather than the committed example, run `riborescue export-web` without
`--sample` over the full `results/amenability_landscape.tsv` and drop the JSON in `public/`.
