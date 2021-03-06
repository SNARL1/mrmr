data {
  int<lower = 1> M;                   // augmented sample size
  int<lower = 1> T;                   // # primary periods
  int<lower = 1> maxJ;                // max # of secondary periods
  int<lower = 0, upper = maxJ> J[T];  // # 2ndary periods for each prim. period
  int<lower = 1> Jtot;                // total number of surveys
  int<lower = 1, upper = T> prim_idx[Jtot];

  // observations
  // 0=NA, 1=not detected, 2=detected
  int<lower = 0, upper = 2> Y[M, T, maxJ];
  int<lower = 0, upper = 1> introduced[M]; // indicator for whether introduced
  int<lower = 0, upper = T> t_intro[M]; // when individuals introduced
  int<lower = 0, upper = 1> removed[M];
  int<lower = 0, upper = T> t_remove[M];

  // index order of surveys (0: NA)
  int<lower = 0, upper = Jtot> j_idx[T, maxJ];
  int<lower = 0, upper = 1> any_surveys[T];

  // fixed effects design matrices
  int<lower = 1> m_detect;
  matrix[Jtot, m_detect] X_detect;
  int<lower = 1> m_surv;
  matrix[M, m_surv] X_surv;

  int<lower = 0, upper = 1> any_recruitment;
}

transformed data {
  int Tm1 = T - 1; // # primary periods - 1
  vector[3] gam_init = [1, 0, 0]';
}

parameters {
  // recruitment
  real alpha_lambda;
  real<lower = 0> sigma_lambda;
  vector[T] eps_lambda;

  // survival
  vector[m_surv] beta_phi;
  real<lower = 0> sigma_phi;
  vector[T] eps_phi;

  // detection params
  vector[m_detect] beta_detect;
}

transformed parameters {
  vector[Jtot] logit_detect;
  vector<lower = 0, upper = 1>[T] lambda;
  vector[M] log_lik;
  matrix[T, 3] forward[M];

  // probability of entering population
  lambda = any_recruitment * inv_logit(alpha_lambda + eps_lambda * sigma_lambda);

  // probability of detection
  logit_detect = X_detect * beta_detect;

  // generate log likelihoods for observation histories
  {
    real acc[3];
    vector[3] gam[T];
    vector[T] phi;
    real ps[3, T, 3];
    real po[3, Jtot, 3];
    real p;
    // s = 1 :: not recruited
    // s = 2 :: alive
    // s = 3 :: dead

    // define probs of state S(t+1) | S(t)
    // first index: S(t)
    // second index: individual
    // third index: t
    // fourth index: S(t + 1)
    for (i in 1:M) {
      // fill in shared values
      for (t in 1:T) {
        phi[t] = inv_logit(X_surv[i, ] * beta_phi + eps_phi[t] * sigma_phi);
        ps[1, t, 3] = 0;       // can't die before being alive
        ps[2, t, 1] = 0;       // can't unenter population
        ps[3, t, 1] = 0;
        ps[2, t, 2] = phi[t];     // survive
        ps[2, t, 3] = 1 - phi[t]; // death
        ps[3, t, 2] = 0; // cannot un-die
        ps[3, t, 3] = 1; // dead stay dead
      }

      if (introduced[i]) {
        // individual has been introduced
        // zero probability of recruiting through t_intro - 1
        for (t in 1:(t_intro[i] - 1)) {
          ps[1, t, 1] = 1;
          ps[1, t, 2] = 0;
          ps[1, t, 3] = 0;
        }

        // timestep of introduction has Pr(recruiting) = 1
        ps[1, t_intro[i], 1] = 0;
        ps[1, t_intro[i], 2] = 1;
        ps[1, t_intro[i], 3] = 0;

        // to avoid NA values, fill in remaining recruitment probs (though they
        // are irrelevant for the likelihood)
        for (t in (t_intro[i] + 1):T) {
          ps[1, t, 1] = 1;
          ps[1, t, 2] = 0;
          ps[1, t, 3] = 0;
        }

      } else {
        for (t in 1:T) {
          ps[1, t, 1] = 1 - lambda[t];
          ps[1, t, 2] = lambda[t];
          ps[1, t, 3] = 0;
        }
      }

      if (removed[i]) {
        if (t_remove[i] < T) {
          ps[2, t_remove[i] + 1, 2] = 0;
          ps[2, t_remove[i] + 1, 3] = 1;
        }
      }

      // observation probabilities
      for (j in 1:Jtot) {
        p = inv_logit(logit_detect[j]);
        po[1, j, 1] = 1;
        po[1, j, 2] = 0;

        if (prim_idx[j] == t_intro[i]) {
          // introductions always happen after surveys, so if an individual is
          // released on primary period t, it has a zero probability of
          // detection
          po[2, j, 1] = 1;
          po[2, j, 2] = 0;
        } else {
          po[2, j, 1] = 1 - p;
          po[2, j, 2] = p;
        }

        po[3, j, 1] = 1;
        po[3, j, 2] = 0;
      }

      for (t in 1:T) {
        for (k in 1:3) {
          for (kk in 1:3) {
            if (t == 1) {
              acc[kk] = gam_init[kk];
            } else {
              acc[kk] = gam[t-1, kk];
            }
            acc[kk] *= ps[kk, t, k];
            if (any_surveys[t]) {
              for (j in 1:J[t]) {
                // only increment the log probability with the likelihood
                // if surveys happened
                acc[kk] *= po[k, j_idx[t, j], Y[i, t, j]];
              }
            }
          }
          gam[t, k] = sum(acc);
        }
        forward[i, t, ] = gam[t, ]';
      }
      log_lik[i] = log(sum(gam[T]));
    } // end loop over individuals
  } // end temporary scope
}

model {
  // priors
  alpha_lambda ~ std_normal();
  sigma_lambda ~ std_normal();
  eps_lambda ~ std_normal();
  beta_detect ~ std_normal();
  beta_phi ~ std_normal();
  sigma_phi ~ std_normal();
  eps_phi ~ std_normal();

  target += sum(log_lik);
}

generated quantities {
  int<lower = 1, upper = 3> s[M, T];  // latent state
  int<lower=0> Nsuper;                // Superpopulation size
  int<lower=0> N[Tm1];                // Actual population size
  int<lower=0> B[Tm1];                // Number of entries

  {
    real ps[3, T, 3];
    vector[T] phi;
    vector[3] tmp;
    // s = 1 :: not recruited
    // s = 2 :: alive
    // s = 3 :: dead

    // define probs of state S(t+1) | S(t)
    // first index: S(t)
    // second index: individual
    // third index: t
    // fourth index: S(t + 1)
    for (i in 1:M) {
      // fill in shared values
      for (t in 1:T) {
        phi[t] = inv_logit(X_surv[i, ] * beta_phi + eps_phi[t] * sigma_phi);
        ps[1, t, 3] = 0;       // can't die before being alive
        ps[2, t, 1] = 0;       // can't unenter population
        ps[3, t, 1] = 0;
        ps[2, t, 2] = phi[t];     // survive
        ps[2, t, 3] = 1 - phi[t]; // death
        ps[3, t, 2] = 0; // cannot un-die
        ps[3, t, 3] = 1; // dead stay dead
      }

      if (introduced[i]) {
        // individual has been introduced
        // zero probability of recruiting through t_intro - 1
        for (t in 1:(t_intro[i] - 1)) {
          ps[1, t, 1] = 1;
          ps[1, t, 2] = 0;
          ps[1, t, 3] = 0;
        }

        // timestep of introduction has Pr(recruiting) = 1
        ps[1, t_intro[i], 1] = 0;
        ps[1, t_intro[i], 2] = 1;
        ps[1, t_intro[i], 3] = 0;

        // to avoid NA values, fill in remaining recruitment probs (though they
        // are irrelevant for the likelihood)
        for (t in (t_intro[i] + 1):T) {
          ps[1, t, 1] = 1;
          ps[1, t, 2] = 0;
          ps[1, t, 3] = 0;
        }

      } else {
        for (t in 1:T) {
          ps[1, t, 1] = 1 - lambda[t];
          ps[1, t, 2] = lambda[t];
          ps[1, t, 3] = 0;
        }
      }

      if (removed[i]) {
        if (t_remove[i] < T) {
          ps[2, t_remove[i] + 1, 2] = 0;
          ps[2, t_remove[i] + 1, 3] = 1;
        }
      }

      // backward sampling
      s[i, T] = categorical_rng(forward[i, T, ]' / sum(forward[i, T]));

      for(t_rev in 1:(Tm1)) {
        int t = T - t_rev;
        int tp1 = t + 1;

        tmp = forward[i, t,]' .* to_vector(ps[, tp1, s[i, tp1]]);
        s[i, t] = categorical_rng(tmp / sum(tmp));
      }
    } // end loop over individuals
  } // end temporary scope

  {
    int al[M, Tm1];
    int d[M, Tm1];
    int alive[M];
    int w[M];

    for (i in 1:M) {
      for (t in 2:T) {
        al[i, t - 1] = (s[i, t] == 2);
      }
      for (t in 1:Tm1)
        d[i, t] = (s[i, t] == al[i, t]);
      alive[i] = sum(al[i]);
    }

    for (t in 1:Tm1) {
      N[t] = sum(al[, t]);
      B[t] = sum(d[, t]);
    }
    for (i in 1:M)
      w[i] = 1 - !alive[i];
    Nsuper = sum(w);
  }

}
