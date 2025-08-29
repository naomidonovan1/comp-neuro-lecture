# comp-neuro-lecture
Slides / materials for computational neuroscience lecture meant for first year neuroscience PhD students.

## Setup
Clone the repo:
```bash
git clone https://github.com/naomidonovan1/comp-neuro-lecture.git
cd comp-neuro-lecture
```

Create a virtual environment (through bash):
```bash
python3.10 -m venv comp-neuro-venv 
source comp-neuro-venv/bin/activate
pip install -r requirements.txt
```

Create a virtual environment (through conda):
```bash
conda create -n comp-neuro-venv python=3.10
conda activate comp-neuro-venv
pip install -r requirements.txt
```

Preview the slides (from your terminal):
First install quarto from https://quarto.org/docs/get-started/
```bash
quarto preview
```