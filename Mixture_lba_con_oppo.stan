functions{
     
     real lba_pdf(real t, real b, real A, real v, real s){
          //PDF of the LBA model
          
          real b_A_tv_ts;
          real b_tv_ts;
          real term_1;
          real term_2;
          real term_3;
          real term_4;
          real pdf;
          
          b_A_tv_ts <- (b - A - t*v)/(t*s);
          b_tv_ts <- (b - t*v)/(t*s);
          term_1 <- v*Phi(b_A_tv_ts);
          term_2 <- s*exp(normal_log(b_A_tv_ts,0,1)); 
          term_3 <- v*Phi(b_tv_ts);
          term_4 <- s*exp(normal_log(b_tv_ts,0,1)); 
          pdf <- (1/A)*(-term_1 + term_2 + term_3 - term_4);
          
          return pdf;
     }
     
     real lba_cdf(real t, real b, real A, real v, real s){
          //CDF of the LBA model
          
          real b_A_tv;
          real b_tv;
          real ts;
          real term_1;
          real term_2;
          real term_3;
          real term_4;
          real cdf;	
          
          b_A_tv <- b - A - t*v;
          b_tv <- b - t*v;
          ts <- t*s;
          term_1 <- b_A_tv/A * Phi(b_A_tv/ts);	
          term_2 <- b_tv/A   * Phi(b_tv/ts);
          term_3 <- ts/A     * exp(normal_log(b_A_tv/ts,0,1)); 
          term_4 <- ts/A     * exp(normal_log(b_tv/ts,0,1)); 
          cdf <- 1 + term_1 - term_2 + term_3 - term_4;
          
          return cdf;
          
     }

     real lba_lpdf(matrix RT, real k, real A, real s, real tau, matrix VALUE, real r, real e){
//hopefully it is lpdf, not lpmf
          
          real t;
          real b;
          real cdf;
          real pdf;

          vector[cols(RT)] prob; //was rows
          real out;
          real prob_neg;
          vector[4] score;
          vector[4] v;

          for (i in 1:cols(RT)){  //was rows

               score <- VALUE[,i];

               for(j in 1:4){
                    v[j] <- r ^ score[j] - 1 + e;
               }
               t <- RT[1,i] - tau;

               if(t > 0){			
                    cdf <- 1;
                    for(j in 1:4){
                         b = A + k;
                         if(RT[2,i] == j){     //
                             pdf <- lba_pdf(t, b, A, v[j], s);
                         }else{
                             cdf <- (1-lba_cdf(t, b, A, v[j], s)) * cdf;
                         }
                    }

                    prob_neg <- 1;
                    for(j in 1:4){
                         prob_neg <- Phi(-v[j]/s) * prob_neg;    
                    }
                    prob[i] <- pdf*cdf;		
                    prob[i] <- prob[i]/(1-prob_neg);	
                    if(prob[i] < 1e-10){
                         prob[i] <- 1e-10;				
                    }
                    
               }else{
                    prob[i] <- 1e-10;			
               }		
          }
          out <- sum(log(prob));
          return out;		
     }

vector lba_lik(matrix RT, real k, real A, real s, real tau, matrix VALUE, real r, real e){
//hopefully it is lpdf, not lpmf

     real t;
     real b;
     real cdf;
     real pdf;

     vector[cols(RT)] prob; //was rows
     vector[cols(RT)] out;
     real prob_neg;
     vector[4] score;
     vector[4] v;

     for (i in 1:cols(RT)){  //was rows

          score <- VALUE[,i];

          for(j in 1:4){
               v[j] <- r ^ score[j] - 1 + e;
          }
          t <- RT[1,i] - tau;

          if(t > 0){
               cdf <- 1;
               for(j in 1:4){
                    b = A + k;
                    if(RT[2,i] == j){     //
                         pdf <- lba_pdf(t, b, A, v[j], s);
                    }else{
                         cdf <- (1-lba_cdf(t, b, A, v[j], s)) * cdf;
                    }
               }

               prob_neg <- 1;
               for(j in 1:4){
                    prob_neg <- Phi(-v[j]/s) * prob_neg;
               }
               prob[i] <- pdf*cdf;
               prob[i] <- prob[i]/(1-prob_neg);
               if(prob[i] < 1e-10){
               prob[i] <- 1e-10;
               }

          }else{
               prob[i] <- 1e-10;
          }
        }
     out <- log(prob);
     return out;
}
     
    matrix lba_rng(real k, real A, real s, real tau, matrix EIG_VALUE, matrix PTS_VALUE, real phi, real r_eig, real r_pts, real e){
          
          int get_pos_drift;	
          int no_pos_drift;
          int get_first_pos;
          vector[4] drift;
          int max_iter;
          int iter;
          real start[4];
          real ttf[4];
          int resp[4];
          real rt;
          matrix[2,20] pred;
          real b;
          vector[4] eig_score;
          vector[4] pts_score;
          vector[4] v;
          vector[cols(EIG_VALUE)] z;

     for (i in 1:cols(EIG_VALUE)){  //was rows
          eig_score <- EIG_VALUE[,i];
          pts_score <- PTS_VALUE[,i];
          z[i] = bernoulli_rng(phi);

          if (z[i] == 1){
              for(j in 1:4){
                   v[j] <- r_eig^eig_score[j]- 1 + e;
              }
          }else{
              for(j in 1:4){
                   v[j] <- r_pts^pts_score[j]- 1 + e;
              }
          }

          get_pos_drift <- 1;
          no_pos_drift <- 0;
          max_iter <- 1000;  // just use for getting a positive drift rate
          iter <- 0;
          while(get_pos_drift){
               for(j in 1:num_elements(v)){
                    drift[j] <- normal_rng(v[j],s);
                    if(drift[j] > 0){
                         get_pos_drift <- 0;
                    }
               }
               iter <- iter + 1;
               if(iter > max_iter){
                    get_pos_drift <- 0;
                    no_pos_drift <- 1;
               }	
          }
          //if both drift rates are <= 0
          //return an infinite response time
          if(no_pos_drift){
               pred[1,i] <- -1;
               pred[2,i] <- -1;
          }else{
               
               for(p in 1:num_elements(v)){
                 	b <- A + k;
                    //start time of each accumulator	
                    start[p] <- uniform_rng(0,A);
                    //finish times
                    ttf[p] <- (b-start[p])/drift[p];
               }
               //rt is the fastest accumulator finish time	
               //if one is negative get the positive drift
               resp <- sort_indices_asc(ttf);
               ttf <- sort_asc(ttf);
               get_first_pos <- 1;
               iter <- 1;
               while(get_first_pos){
                    if(ttf[iter] > 0){
                         pred[1,i] <- ttf[iter] + tau;
                         pred[2,i] <- resp[iter];
                         get_first_pos <- 0;
                    }
                    iter <- iter + 1;
               }
               }
        }
          return pred;	
     }

}

data {
     int NUM_NODES;
     matrix[4,20] EIG_VALUE_e;
     matrix[4,20] PTS_VALUE_e;
     matrix[4,20] EIG_VALUE_n;
     matrix[4,20] PTS_VALUE_n;
     int NUM_SUBJ;
     matrix[2,20] RT[NUM_SUBJ*2];
}
//transformed data{
  //   matrix[2,20] RT_e[NUM_SUBJ];
   //  matrix[2,20] RT_n[NUM_SUBJ];
 //    RT_e = RT[1:NUM_SUBJ];
   //  RT_n = RT[NUM_SUBJ+1:NUM_SUBJ*2];
//}

parameters {
     real<lower=0> k_e[NUM_SUBJ];
     real<lower=0> k_n[NUM_SUBJ];
     real<lower=0> A[NUM_SUBJ];
     real<lower=0> tau[NUM_SUBJ];

     real<lower=0.01> r_eig_e[NUM_SUBJ];
     real<lower=0.01> r_pts_e[NUM_SUBJ];
     real<lower=0.01> r_eig_n[NUM_SUBJ];
     real<lower=0.01> r_pts_n[NUM_SUBJ];

     real<lower=0> e[NUM_SUBJ]; // single e

     real<lower=0,upper=1> phi_e[NUM_SUBJ];
     real<lower=0,upper=1> phi_n[NUM_SUBJ];

     real<lower=0> k_e_sigma;
     real<lower=0> k_n_sigma;
     real<lower=0> A_sigma;

     real<lower=0> k_e_mu;
     real<lower=0> k_n_mu;
     real<lower=0> A_mu;

     real<lower=0> r_eig_e_sigma;
     real<lower=0> r_eig_e_mu;
     real<lower=0> r_pts_e_sigma;
     real<lower=0> r_pts_e_mu;

     real<lower=0> r_eig_n_sigma;
     real<lower=0> r_eig_n_mu;
     real<lower=0> r_pts_n_sigma;
     real<lower=0> r_pts_n_mu;

     real<lower=0> e_sigma;
     real<lower=0> e_mu;

     real<lower=0> phi_e_k;
     real<lower=0,upper=1> phi_e_mu;
     real<lower=0> phi_n_k;
     real<lower=0,upper=1> phi_n_mu;

}

transformed parameters {
     real s;
     real<lower=0> phi_e_a;
     real<lower=0> phi_e_b;
     real<lower=0> phi_n_a;
     real<lower=0> phi_n_b;

     s <- 1;
     phi_e_a <- phi_e_k * phi_e_mu;
     phi_e_b <- phi_e_k * (1-phi_e_mu);
     phi_n_a <- phi_n_k * phi_n_mu;
     phi_n_b <- phi_n_k * (1-phi_n_mu);
}

model {

     A_mu ~ uniform(.01,50)T[0,];
     k_e_mu ~ uniform(.01,20)T[0,];
     k_n_mu ~ uniform(.01,20)T[0,];

     A_sigma ~ gamma(1,1);
     k_e_sigma ~ gamma(1,1);
     k_n_sigma ~ gamma(1,1);

     r_eig_e_mu ~ uniform(.01,30)T[0,];
     r_eig_e_sigma ~ gamma(1,1);
     r_pts_e_mu ~ uniform(.01,30)T[0,];
     r_pts_e_sigma ~ gamma(1,1);

     r_eig_n_mu ~ uniform(.01,30)T[0,];
     r_eig_n_sigma ~ gamma(1,1);
     r_pts_n_mu ~ uniform(.01,30)T[0,];
     r_pts_n_sigma ~ gamma(1,1);

     e_mu ~ uniform(.01,20)T[0,];
     e_sigma ~ gamma(1,1);

     phi_e_k ~ gamma(10,1);
     phi_e_mu ~ beta(1,1);

     phi_n_k ~ gamma(10,1);
     phi_n_mu ~ beta(1,1);

     for(i in 1:NUM_SUBJ){
          A[i] ~ normal(A_mu,A_sigma)T[0,];
          tau[i] ~ uniform(.01,10)T[0,];
          k_e[i] ~ normal(k_e_mu,k_e_sigma)T[0,];
          k_n[i] ~ normal(k_n_mu,k_n_sigma)T[0,];

          r_eig_e[i] ~ normal(r_eig_e_mu,r_eig_e_sigma)T[0.01,];
          r_pts_e[i] ~ normal(r_pts_e_mu,r_pts_e_sigma)T[0.01,];
          r_eig_n[i] ~ normal(r_eig_n_mu,r_eig_n_sigma)T[0.01,];
          r_pts_n[i] ~ normal(r_pts_n_mu,r_pts_n_sigma)T[0.01,];
          e[i] ~ normal(e_mu,e_sigma)T[0,];

          phi_e[i] ~ beta(phi_e_a,phi_e_b);
          phi_n[i] ~ beta(phi_n_a,phi_n_b);
          //phi_e[i] ~ uniform(0,1);
          //phi_n[i] ~ uniform(0,1);
          target += log_mix(phi_e[i],
                            lba_lpdf(RT[i] | k_e[i],A[i],s,tau[i],EIG_VALUE_e,r_eig_e[i],e[i]),
                            lba_lpdf(RT[i] | k_e[i],A[i],s,tau[i],PTS_VALUE_e,r_pts_e[i],e[i]));
          target += log_mix(phi_n[i],
                            lba_lpdf(RT[NUM_SUBJ+i] | k_n[i],A[i],s,tau[i],EIG_VALUE_n,r_eig_n[i],e[i]),
                            lba_lpdf(RT[NUM_SUBJ+i] | k_n[i],A[i],s,tau[i],PTS_VALUE_n,r_pts_n[i],e[i]));

     }

}


generated quantities {
     matrix[2,20] pred[NUM_SUBJ*2];
     vector[20] log_eig[NUM_SUBJ*2];
     vector[20] log_pts[NUM_SUBJ*2];
     vector[20] log_lik[NUM_SUBJ*2];

     for (i in 1:NUM_SUBJ){
          pred[i] <- lba_rng(k_e[i],A[i],s,tau[i],EIG_VALUE_e,PTS_VALUE_e,phi_e[i],r_eig_e[i],r_pts_e[i],e[i]);
          pred[NUM_SUBJ+i] <- lba_rng(k_n[i],A[i],s,tau[i],EIG_VALUE_n,PTS_VALUE_n,phi_n[i],r_eig_n[i],r_pts_n[i],e[i]);

          log_eig[i] = lba_lik(RT[i],k_e[i],A[i],s,tau[i],EIG_VALUE_e,r_eig_e[i],e[i]);
          log_pts[i] = lba_lik(RT[i],k_e[i],A[i],s,tau[i],PTS_VALUE_e,r_pts_e[i],e[i]);
          log_eig[NUM_SUBJ+i] = lba_lik(RT[NUM_SUBJ+i],k_n[i],A[i],s,tau[i],EIG_VALUE_n,r_eig_n[i],e[i]);
          log_pts[NUM_SUBJ+i] = lba_lik(RT[NUM_SUBJ+i],k_n[i],A[i],s,tau[i],PTS_VALUE_n,r_pts_n[i],e[i]);
          for (j in 1:20){
               log_lik[i,j] <- log_sum_exp(log(phi_e[i]) + log_eig[i,j],
                                           log(1-phi_e[i]) + log_pts[i,j]);
               log_lik[NUM_SUBJ+i,j] <- log_sum_exp(log(phi_n[i]) + log_eig[NUM_SUBJ+i,j],
                                                    log(1-phi_n[i]) + log_pts[NUM_SUBJ+i,j]);
          }
     }
}
