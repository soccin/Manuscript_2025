project_name/
├── dat/
│   ├── raw/
│   ├── proc/
│   └── meta/
├── src/
├── out/
│   ├── fig/          # generated figures
│   ├── tab/          # generated tables
│   └── qc/           # quality control reports
├── doc/
│   ├── methods/      # detailed methods documentation
│   ├── notes/        # lab notebook style notes
│   └── refs/         # bibliography, citations
├── slides/           # or pres/ for presentations
│   ├── lab_meet/
│   ├── conf/
│   └── collab/
├── ms/               # manuscript
│   ├── draft/
│   ├── fig/          # manuscript figures (symlinks to out/fig/)
│   ├── supp/         # supplementary materials
│   └── rev/          # reviews and responses
└── README
```

## File Organization Strategy

**doc/methods/**
- One file per major analysis step
- Clear links to code: `rnaseq_align.md` → references `src/align.sh`
- Include command lines, parameters, software versions

**slides/**
- Organize by audience/venue
- Keep source files (`.pptx`, `.key`, or better yet `.md` for reveal.js)

**ms/**
- `draft/` - working manuscript files
- `fig/` - polished figures for manuscript (often symlinks: `ln -s ../../out/fig/fig1.pdf fig/fig1.pdf`)
- `supp/` - supplementary tables, figures, methods
- `rev/` - reviewer comments and response letters (when you get there)

## Alternative: More Concise

If you want even shorter top-level directories:
```
project_name/
├── dat/
├── src/
├── out/
├── doc/
├── talk/        # instead of slides/
└── paper/       # instead of ms/
