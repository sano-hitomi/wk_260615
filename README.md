# Project Name

> **This is a template repository.** Use this as a starting point for new projects. Replace this line with a 1–2 sentence overview of your project.

---

## Directory Structure

```
.
├── data/
│   ├── dummy/              # Dummy data (tracked by Git)
│   │   ├── raw/            #   Raw data before processing
│   │   └── processed/      #   Cleaned and processed data
│   └── production/         # Production data (excluded from Git via .gitignore)
│       ├── raw/            #   Raw data before processing
│       └── processed/      #   Cleaned and processed data
├── src/
│   ├── python/             # Python scripts
│   └── r/                  # R scripts
├── notebooks/              # Exploratory analysis (Jupyter Notebook / R Markdown)
├── output/
│   ├── figures/            # Figures and plots
│   └── reports/            # Reports and summary results
├── docs/                   # Documentation and specifications
├── config/                 # Configuration files (do NOT include secrets or API keys)
├── .gitignore
└── README.md
```

---

## Data Handling Rules

| Folder | Git Tracked | Purpose |
|---|---|---|
| `data/dummy/` | ✅ Yes | Dummy data for development and testing |
| `data/production/` | ❌ No | Actual production data (do not share or publish) |

> **Note:** `data/production/` is excluded via `.gitignore`.
> Be careful not to accidentally commit production data.

---

## Git Setup (for new projects based on this template)

```bash
git init
git config core.hooksPath .githooks
git add .
git commit -m "Initial commit: project structure"
```

> **Note:** `git config core.hooksPath .githooks` enables the pre-commit hook that prevents files larger than 50MB from being committed.

---

## Contact

- Maintainer:
- Email:
