# Magma Ocean–Atmosphere Model

The main script is:

```text
convection_redox.py
```

This program models the structure of an atmosphere overlying a magma ocean for a range of magma ocean potential temperatures. It calculates the atmospheric thermal structure and composition after volatile degassing from the magma ocean.

## Setting the Initial Composition

The initial volatile inventory is defined in **`convection_redox.py`** (around lines 120–123).

### 1. Initial Carbon and Water Contents

```python
# IW+4
xC_base, xH_base = 2.e-4, 4.e-4
```

These values specify the initial concentrations of carbon and water in the magma ocean before degassing.

- `xC_base`: total carbon concentration
- `xH_base`: total hydrogen (water) concentration

### 2. Carbon Speciation

```python
xCO_init, xCO2_init = xC_base * 28/44 * 0., xC_base * 1.
```

These lines determine how the total carbon is partitioned between **CO₂** and **CO**.

In the example above:

- 100% of the carbon is initially assigned to **CO₂**.
- 0% of the carbon is initially assigned to **CO**.

### 3. Hydrogen Speciation

```python
xH2_init, xH2O_init = xH_base * 2/18 * 0., xH_base * 1.
```

These lines determine how the total hydrogen is partitioned between **H₂O** and **H₂**.

In the example above:

- 100% of the hydrogen is initially assigned to **H₂O**.
- 0% of the hydrogen is initially assigned to **H₂**.

## Example Configuration

The settings shown above correspond to an **IW+4** redox state with:

- Carbon concentration: `2 × 10⁻⁴`
- Water concentration: `4 × 10⁻⁴`
- All carbon initially as **CO₂**
- All hydrogen initially as **H₂O**
