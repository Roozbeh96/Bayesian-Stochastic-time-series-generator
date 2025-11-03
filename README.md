# Bayesian-Stochastic generation of (spatio/temporal)-series signal

Tool to stochastically generate (spatio/temporal)-series velocity signal. This velocity signal can be used in stochastic generation of turbulent boundary layer velocity field.

## Contents
- [python package](#python-package)
- [Algorithm](#algorithm)

## Python package
<!-- brief blurb or link to docs/install/usage -->
<details>
<summary> Installation:</summary>

1. 🛠️ Installing Poetry

To install [Poetry](https://python-poetry.org/) (Python dependency management and packaging tool), run the following command in your terminal:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

After installation, make sure Poetry is in your `PATH`
- macOS/Linux
```bash
export PATH="$HOME/.local/bin:$PATH
```

Verify installation:

```bash
poetry version
```

Keep venv inside the project (works great with VS Code) poetry 

```bash
config virtualenvs.in-project true
```

2. Installing environment: 

This environment is set with Python 3.13. Change the requires-python = ">=3.13" in pyproject.tmol file if you have other versions on your PC.

run:
```bash
poetry install
```

In case if you want to make environment from scratch, run:(Do not recommended)

```bash
poetry new project_name
```
If you want to delete the caches use:
```bash
git rm -r --cached .
```

Then re download the .git file and do:
```bash
git add .
git commit -m 'init'
git push
```
</details>

<details>
  <summary>Click to see the documentation for Python files</summary>

### Python files
`Bayesian_ST/src/bayesian_st`<br>
-`main.py:` Define the parameters.<br>


`Bayesian_ST/src/bayesian_st/utils`<br>
-`Stochastic_generation.py:` Defining object.<br>
-`Short_time_analysis.m:` Chope signal into the desired length. For this file, you need experimental dataset, which I have used Hotwire time resolved dataset ([Link])(https://conservancy.umn.edu/items/e2f507c9-570d-46b6-b70c-939877caf668). Fig1 and Fig2 can be generated using this file. There is no need to run this file, just use the output vectors, which are standard deviation of fluctuating velocity signal (u^{\prime}/u_{\tau}, w^{\prime}/u_{\tau}) of the chopped time-series velocity signal (,std_uprime_chopped.mat,std_wprime_chopped.mat). At the end of this file, we only use the results of the standard deviation of the chopped HotWT7. To justify our claim, I did hypothesis testing.<br>
$$
H_{0}: \abs(\sigma_{\text{chopped_WT7}}-\sigma_{\text{chopped_WT10}}) < \delta 
H_{1}: \abs(\sigma_{\text{chopped_WT7}}-\sigma_{\text{chopped_WT10}}) > \delta
$$





</details>

## Algorithm
Fig 1: The velocity signal is chopped into the desired length scale.<br>
<img src="Fig/Fig17.png" alt="what!!!" width="400"/><br>
Fig 2: Variability of the velocity signal after chopping for different datasets.<br>
<img src="Fig/Fig2.png" alt="what!!!" width="400"/><br>
Fig 3: Shows the algorithm of generating velocity signal.<br>
<img src="Fig/Fig20.png" alt="what!!!" width="400"/><br>
<img src="Fig/Fig24.png" alt="what!!!" width="400"/><br>
<img src="Fig/Fig27.png" alt="what!!!" width="400"/><br>
<img src="Fig/Fig28.png" alt="what!!!" width="400"/><br>