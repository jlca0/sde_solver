# SDE: Euler-Maruyama y Milstein en C++ con pybind11

## Compilar
```bash
pip install -r requirements.txt
cmake -S . -B build -Dpybind11_DIR=$(python -m pybind11 --cmakedir)
cmake --build build
```

## Ejecutar
```bash
cd py
python compare.py
```