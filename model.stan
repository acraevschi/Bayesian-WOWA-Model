functions{


    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }


    matrix cov_GPL1(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data{
    int num_doc;
    int pos[11272];
    int weight[11272];
    int textID[11272];
    int flag[11272];
    int role[11272];
    int doculect[11272];
    matrix[24,24] phylo;
    matrix[24,24] geo_dist;
}
parameters{
    matrix[2,5] z_rf;
    vector<lower=0>[2] sigma_rf;
    cholesky_factor_corr[2] L_Rho_rf;
    vector[196] text;
    real beta_w;
    vector[num_doc] phy_z;
    real<lower=0> etasq_phy;
    real<lower=0> rhosq_phy;
    real<lower=0> delta_phy;
    vector[num_doc] geo_z;
    real<lower=0> etasq_geo;
    real<lower=0> rhosq_geo;
    real<lower=0> delta_geo;
}
transformed parameters{
    matrix[5,2] role_flag;
    vector[num_doc] phy;
    matrix[num_doc,num_doc] phy_Sigma;
    matrix[num_doc,num_doc] phy_SIGMA;
    vector[num_doc] geo;
    matrix[num_doc,num_doc] geo_Sigma;
    matrix[num_doc,num_doc] geo_SIGMA;
    geo_SIGMA = cov_GPL2(geo_dist, etasq_geo, rhosq_geo, delta_geo);
    geo_Sigma = cholesky_decompose(geo_SIGMA);
    geo = geo_Sigma * geo_z;
    phy_SIGMA = cov_GPL1(phylo, etasq_phy, rhosq_phy, delta_phy);
    phy_Sigma = cholesky_decompose(phy_SIGMA);
    phy = phy_Sigma * phy_z;
    role_flag = (diag_pre_multiply(sigma_rf, L_Rho_rf) * z_rf)';
}
model{
    vector[11272] p;
    delta_geo ~ exponential( 2 );
    rhosq_geo ~ exponential( 0.35 );
    etasq_geo ~ exponential( 0.35 );
    geo_z ~ normal( 0 , 1 );
    delta_phy ~ exponential( 2 );
    rhosq_phy ~ exponential( 0.35 );
    etasq_phy ~ exponential( 0.35 );
    phy_z ~ normal( 0 , 1 );
    beta_w ~ normal( 0 , 1.5 );
    text ~ normal( -2 , 0.5 );
    L_Rho_rf ~ lkj_corr_cholesky( 2 );
    sigma_rf ~ exponential( 1 );
    to_vector( z_rf ) ~ normal( 0 , 1 );
    for ( i in 1:11272 ) {
        p[i] = phy[doculect[i]] + geo[doculect[i]] + role_flag[role[i], flag[i]] + text[textID[i]] + beta_w * weight[i];
        p[i] = inv_logit(p[i]);
    }
    pos ~ binomial( 1 , p );
}
generated quantities{
    vector[11272] log_lik;
    vector[11272] p;
    matrix[2,2] Rho_rf;
    Rho_rf = multiply_lower_tri_self_transpose(L_Rho_rf);
    for ( i in 1:11272 ) {
        p[i] = phy[doculect[i]] + geo[doculect[i]] + role_flag[role[i], flag[i]] + text[textID[i]] + beta_w * weight[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:11272 ) log_lik[i] = binomial_lpmf( pos[i] | 1 , p[i] );
}

