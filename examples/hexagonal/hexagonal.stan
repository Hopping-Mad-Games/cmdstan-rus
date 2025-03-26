// Hexagonal single crystal elastic constants estimation model
functions {
  vector mech_init(int P, real X, real Y, real Z, real density);
  matrix mech_rotate(matrix C, vector q);
  vector mech_rus(int N, vector lookup, matrix C);
  vector cu2qu(vector cu);
}

// Input data
data {
  int<lower = 1> P; // Order of polynomials for Rayleigh-Ritz approx
  int<lower = 1> N; // Number of resonance modes

  // Sample dimensions and physical properties
  real<lower = 0.0> X; // width
  real<lower = 0.0> Y; // length
  real<lower = 0.0> Z; // height
  real<lower = 0.0> density;

  // Measured resonance frequencies
  vector[N] y;
}

transformed data {
  int L = (P + 1) * (P + 2) * (P + 3) / 6;
  vector[L * L * 3 * 3 * 21] lookup;

  lookup = mech_init(P, X, Y, Z, density);
}

// Parameters to estimate - Hexagonal crystals have 5 independent elastic constants
parameters {
  real<lower = 0.0> c11;  // in-plane stiffness
  real<lower = 0.0> c33;  // out-of-plane stiffness
  real<lower = 0.0> c44;  // out-of-plane shear modulus
  real<lower = 0.0> c12;  // in-plane Poisson ratio effect
  real<lower = 0.0> c13;  // coupling between in-plane and out-of-plane deformations
  real<lower = 0.0> sigma; // measurement noise
  
  // Crystal orientation parameters (optionally include if rotations are needed)
  vector<lower = -1.0725146985555127, upper = 1.0725146985555127>[3] cu; // rotation between sample & crystal axes
}

// Build a 6x6 stiffness matrix for hexagonal symmetry and apply rotation
transformed parameters {
  matrix[6, 6] C;
  vector[4] q;
  real c66 = (c11 - c12) / 2.0; // dependent elastic constant
  real zero = 0.0;
  
  // Initialize matrix to zeros
  for (i in 1:6)
    for (j in 1:6)
      C[i, j] = zero;
  
  // Fill in the elastic constants for hexagonal symmetry
  C[1, 1] = c11;
  C[2, 2] = c11;
  C[3, 3] = c33;
  C[4, 4] = c44;
  C[5, 5] = c44;
  C[6, 6] = c66;
  
  C[1, 2] = c12;
  C[2, 1] = c12;
  
  C[1, 3] = c13;
  C[3, 1] = c13;
  
  C[2, 3] = c13;
  C[3, 2] = c13;

  // Convert Cayley-Klein parameters to quaternions and rotate stiffness matrix
  q = cu2qu(cu);
  C = mech_rotate(C, q);
}

// Probabilistic model
model {
  // Priors on elastic constants - adjust these based on expected ranges for your material
  c11 ~ normal(100.0, 50.0);
  c33 ~ normal(100.0, 50.0);
  c44 ~ normal(50.0, 25.0);
  c12 ~ normal(50.0, 25.0);
  c13 ~ normal(40.0, 20.0);
  
  // Prior on noise level (in kHz)
  sigma ~ normal(0, 0.2);
  
  // Likelihood: measured frequencies are normally distributed around predicted frequencies
  y ~ normal(mech_rus(N, lookup, C), sigma);
}

generated quantities {
  vector[N] yhat;     // posterior predictive samples
  vector[N] errors;   // residuals
  real anisotropy_ratio = c33 / c11;  // c-axis vs a-axis stiffness ratio

  {
    vector[N] freqs = mech_rus(N, lookup, C);
    for(n in 1:N) {
      errors[n] = y[n] - freqs[n];
      yhat[n] = normal_rng(freqs[n], sigma);
    }
  }
}
