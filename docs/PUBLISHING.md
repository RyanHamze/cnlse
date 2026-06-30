# Publishing The Public GitHub Repository

The repository is prepared locally as `cnlse` with `main` as the default branch.

## Preferred One-Time Publish

Install and authenticate the GitHub CLI, then run these commands from the repository root:

```bash
gh auth login
gh repo create RyanHamze/cnlse --public --source=. --remote=origin --push --description "Fiber optics modulation simulations using NLSE/CNLSE and split-step Fourier methods."
```

## If The Empty GitHub Repo Already Exists

If `RyanHamze/cnlse` has already been created on GitHub, run:

```bash
git remote add origin git@github.com:RyanHamze/cnlse.git
git push -u origin main
```

If SSH is not configured, use HTTPS instead:

```bash
git remote add origin https://github.com/RyanHamze/cnlse.git
git push -u origin main
```

## Current Local Verification

Before publishing, this repo was verified with:

```bash
python -m pytest tests -q
python -m py_compile src/cnlse/*.py tests/*.py examples/ook_100km_compensated.py app.py
python examples/ook_100km_compensated.py
python -m pip install -e .
cnlse-example
```

## Repository Metadata

- Name: `cnlse`
- Title: `Fiber Optics Modulation Using Coupled Nonlinear Schrodinger Equation`
- License: MIT
- Primary author/maintainer: Ryan Maher Hamze
- Original 2009 MATLAB project authors: Ryan Maher Hamze, Krayem Ahmad Ajjawi, Meriam Hawwa
- Related ebook ASIN: `B077Z6FMF3`

## Streamlit Community Cloud

Deploy from:

- Repository: `RyanHamze/cnlse`
- Branch: `main`
- Main file path: `streamlit_app.py`

The root `requirements.txt` installs the local package plus Streamlit app dependencies.
