# Hexagonal Single Crystal Models

This directory contains Stan models and example data files for estimating the elastic constants of hexagonal single crystals using Resonant Ultrasound Spectroscopy (RUS).

## Hexagonal Crystal Symmetry

Hexagonal crystals have 5 independent elastic constants:
- `c11`: In-plane stiffness (a-axis direction)
- `c33`: Out-of-plane stiffness (c-axis direction)
- `c44`: Out-of-plane shear modulus
- `c12`: In-plane Poisson ratio effect
- `c13`: Coupling between in-plane and out-of-plane deformations

The sixth constant `c66` is related to these by: `c66 = (c11 - c12)/2`

The full stiffness matrix for a hexagonal crystal has the form:
```
[c11 c12 c13  0   0   0  ]
[c12 c11 c13  0   0   0  ]
[c13 c13 c33  0   0   0  ]
[ 0   0   0  c44  0   0  ]
[ 0   0   0   0  c44  0  ]
[ 0   0   0   0   0  c66 ]
```

## Models Provided

Two models are provided:

1. `hexagonal.stan`: Full model that includes rotation parameters to account for misalignment between the sample axes and crystal axes. This model is more realistic as most samples aren't perfectly aligned with the crystal axes.

2. `hexagonal_simple.stan`: Simplified model that assumes the sample axes are aligned with the crystal axes. Use this as a starting point if you're sure your sample is well-aligned or if you want a faster, simpler model.

## Data Files

Example data files are provided in two formats:

- `hexagonal.data.json`: JSON format for use with CmdStan directly 
- `hexagonal.dat`: R-style data format

The data represents a hypothetical titanium sample with density around 4500 kg/m³ and dimensions of approximately 6.5 × 7.4 × 11.2 mm.

## Example Usage

To run the model with CmdStan:

```bash
./path/to/cmdstan/bin/stanc examples/hexagonal/hexagonal.stan -o examples/hexagonal/hexagonal.hpp
make examples/hexagonal/hexagonal
./examples/hexagonal/hexagonal sample data file=examples/hexagonal/hexagonal.data.json output file=output.csv
```

For the simple model:

```bash
./path/to/cmdstan/bin/stanc examples/hexagonal/hexagonal_simple.stan -o examples/hexagonal/hexagonal_simple.hpp
make examples/hexagonal/hexagonal_simple
./examples/hexagonal/hexagonal_simple sample data file=examples/hexagonal/hexagonal.data.json output file=output_simple.csv
```

## Interpreting Results

The model outputs:
- Estimates of the 5 independent elastic constants (`c11`, `c33`, `c44`, `c12`, `c13`)
- The measurement noise parameter (`sigma`)
- For the rotated model, the orientation parameters (`cu`)
- Anisotropy ratio (`c33/c11`): Values > 1 indicate the crystal is stiffer along the c-axis than in the basal plane

## Common Hexagonal Materials

Typical room temperature elastic constants for common hexagonal materials (in 10¹¹ Pa):

| Material | c11 | c33 | c44 | c12 | c13 |
|----------|-----|-----|-----|-----|-----|
| Ti       | 1.62| 1.81| 0.47| 0.92| 0.69|
| Zn       | 1.65| 0.62| 0.39| 0.31| 0.50|
| Zr       | 1.43| 1.65| 0.32| 0.73| 0.65|
| Mg       | 0.59| 0.62| 0.16| 0.26| 0.22|
| Be       | 2.92| 3.36| 1.63| 0.27| 0.14|

The priors in the models are set for titanium-like materials but can be adjusted for other materials as needed.
