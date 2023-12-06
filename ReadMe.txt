The steps below will allow one to implement the policy search algorithm, as described in the paper: Guan Q, Reich BJ, Laber E and Bandyopadhyay D. (2020). Bayesian nonparametric policy search with applications to periodontal recall intervals, Journal of the American Statistical Association - Applications & Case Studies,  105(531), 1066-1078

1. "Sim_DPM" is the main file used to run the demo code. It consists of five main steps:
   (1) Set simulation parameters
   (2) Generate a dataset
   (3) Fit the model
   (4) Policy search
   (5) Store the output

2. "SimDesign" sets parameters for simulating one dataset, including coefficients and variance in the model, group assignment probablility, and sample size.

3. "Sim_One_Data" simulates one dataset. It simulates the disease progreesion within 60 months for each subject, and computes the sufficient statistics that will be used in MCMC sampling.

4. "MCMC_DPM" contains the function of fitting the DPM model. The function takes the sufficient statistics, number of mixture components, hyperparameters, 
    number of iterations in the MCMC sampling and number of burn-in samples as the input, and the posterior samples of the parameters as the output.

5. "optimization" implements the policy search algorithm to optimize the feature weights that maximize the value.

6. "function" defines functions that are needed in above procedures.
   (1) "gen.one.sub": simulates the disease progression within 60 months for one subject.
	Arguments:
	  params: parameters of the model;
	  policy: the policy used to recommend the recall interval. It can be "riskscore", "randompolicy" or "policy6";
	  alpha: the features weights in the risk score policy.
	Value produced: 
	  a list object containing simulated disease progression for one subject.

   (2) "gen.m.subs": simulates the disease progression within 60 months for m subjects.
	Arguments:
	  m: the number of subjects to simulate;
	  params: model parameters;
	  policy: the policy used to recommend the recall interval, it can be "riskscore", "randompolicy" or "policy6";
	  alpha: the features weights in the risk score policy;
	  utility: utility function, it can be "red" or "ave";
          indX_policy: indicator of the elements in X that will be a feature in the policy.
	Value produced: 
	  a list object containing simulated disease progression for m subjects.

   (3) "get.sufficient.stats": computes the sufficient statistics X'X, X'Y, Y'Y.
	Arguments:
	  X: baseline covariates;
	  Y: disease response trajectory;
	  delta: actual time between visits;
          X1: all independent variables in the compliance model;
	  X2: all independent variables in the disease progression model.
	Value produced: 
	  a list object containing the sufficient statistics.

   (4) "randompolicy": policy that recommends recall interval A[t] based on logit(Pr(A[t]=3))=Y[t-1].
	Arguments:
	  X: baseline covariates that are used as a feature;
	  DA: compliance;
          Y: current disease status;
	  diffY: disease progression since last visit;
	  alpha: feature weights.
	Value produced: 
	  recommended recall interval.

   (5) "policy6": baseline policy that recommends A=6 months between visits for all subjects and t.
	Arguments:
	  X: baseline covariates that are used as a feature;
	  DA: compliance;
          Y: current disease status;
	  diffY: disease progression since last visit;
	  alpha: feature weights;
	Value produced: 
	  recommended recall interval.

   (6) "riskscore": policy based on risk score.
	Arguments:
	  X: baseline covariates that are used as a feature;
	  DA: compliance;
          Y: current disease status;
	  diffY: disease progression since last visit;
	  alpha: feature weights;
	Value produced: 
	  recommended recall interval.

   (7) "get.value": computes the value of a policy defined by a.
	Arguments:
	  a: feature weights;
	  params: model parameters;
          m: number of subjects;
	  policy: the policy used to recommend the recall interval, it can be "riskscore", "randompolicy" or "policy6";
	  thresh: threshold used in the risk score policy;
	  utility: utility function, it can be "red" or "ave";
          indX_policy: indicator of the elements in X that will be a feature in the policy. 
	Value produced: 
	  the value and average recall time of a policy.


   (8) "get.thresh": estimates the threshold satisfying the average recall time constraint for a policy defined by a.
	Arguments:
	  a: feature weights;
	  params: model parameters;
          m: number of subjects;
	  policy: the policy used to recommend the recall interval, it can be "riskscore", "randompolicy" or "policy6";
	  npts: the number of candidate points to search over in each round;
          indX_policy: indicator of the elements in X that will be a feature in the policy. 
	Value produced: 
	  the threshold satisfying the average recall time constraint.

   (9) "ccd_sphere": central composite design on a sphere.
	  width: the width of the design;
	  p: number of dimensions;
	  m: number of grid points for one dimension;
	Value produced: 
	  a matrix of central composite design.

   (10) "GP_MLE": gets maximum likelihood estimates of parameters in the Gaussian process
	  Y: a vector of observed responses;
	  X: a matirx of observed independent variables.
	Value produced: 
	  maximum likelihood estimates of parameters in the Gaussian process.

   (11) "GP_predict": predicts the values at new data points in the Gaussian process.
	  Xnew: a matrix of new independent variables;
	  Y: a vector of observed responses;
	  X: a matirx of observed independent variables;
	  parmas: Gaussian process parameters.
	Value produced: 
	  a list project containing the predicted mean and standard deviation of values at new data points in the Gaussian process.

   (12) "gain": computes the expected decrease.
	  thresh: current minimum;
	  mu: predicted mean;
	  sigma: predicted standard deviation.
	Value produced: 
	  the expected decrease.

7. "binary_covariate" folder contains all the codes and functions for the simulation study with binary covariates.
   It includes similar functions and steps, but has a more general version of "MCMC_DPM" function which can deal with binary covariates.

8. "model_misspecified" folder contains all the codes and functions for the simulation study for model misspecification.
   It includes functions that estimate the value associated with a policy using the true generative model and the proposed model respectively.
   Other functions and steps are similar as that in the main simulation study.










