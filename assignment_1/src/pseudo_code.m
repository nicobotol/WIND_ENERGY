load airfoid data
load blade parameters
initialize a vector of lambda
initialize a vector of Theta_p

% Main code
loop over lambda
  loop over Theta_p
    loop over blade radius
      run BEM code
    stop looping over r
    integrate the partial results of different r
    store lambda, Theta_p, cP, cT
  stop looping over Theta_p
stop looping over lambda 

% BEM code
initialize a, a'
  loop up to a max number of iterations
    update a, a'
    compute Phi
    compute Theta
    compute alpha
    interpolate cl and cd from the table
    compute cn, ct
    compute F
    compute a better guess for a, a'
    if |a - a_old| < epsilon & |a' - a'_old| < epsilon
      break
    else
      go on
    end
  stop looping



  initialize a velocity vector
  initilize theta_p vector
  loop over velocities
    compute omega
    compute lambda
    loop over pitch angle
      run BEM code
      compute P(theta)
    stop looping over pitch angles
    find theta_p s.t P(theta_p) = Prated
    store V0 and theta_p
  stop looping over velocities
