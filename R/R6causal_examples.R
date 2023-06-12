#' @include R6causal.R 
NULL

#' SCM "backdoor" used in the examples.
#'
#' Variable z fulfills the back-door criterion for P(y|do(x))
#' @examples
#' backdoor
#' backdoor$plot()
#' @export
backdoor <- SCM$new("backdoor",
                    uflist = list(
                      uz = function(n) {return(stats::runif(n))},
                      ux = function(n) {return(stats::runif(n))},
                      uy = function(n) {return(stats::runif(n))}
                    ),
                    vflist = list(
                      z = function(uz) {
                        return(as.numeric(uz < 0.4))},
                      x = function(ux, z) {
                        return(as.numeric(ux < 0.2 + 0.5*z))},
                      y = function(uy, z, x) {
                        return(as.numeric(uy < 0.1 + 0.4*z + 0.4*x))}
                    )
)

#' SCM "frontdoor" used in the examples.
#'
#' Variable z fulfills the front-door criterion for P(y|do(x))
#' @examples
#' frontdoor
#' frontdoor$plot()
#' @export
frontdoor <- SCM$new("frontdoor",
                     uflist = list(
                       uz = function(n) {return(stats::runif(n))},
                       ux = function(n) {return(stats::runif(n))},
                       uy = function(n) {return(stats::runif(n))},
                       u = function(n) {return(stats::runif(n))}
                     ),
                     vflist = list(
                       x = function(ux,u) {
                         return( as.numeric(ux < 0.2 + 0.7*u))},
                       z = function(uz,x) {
                         return( as.numeric(uz < 0.1 + 0.7*x))},
                       y = function(uy,z,u) {
                         return( as.numeric(uy < 0.1 + 0.3*z + 0.5*u + 0.2*u*z))}
                     )
)

#' SCM "trapdoor" used in the examples.
#'
#' Variable z is a trapdoor variable for P(y|do(x))
#' @references
#' J. Helske, S. Tikka, J. Karvanen (2021). Estimation of causal effects with
#' small data in the presence of trapdoor variables,
#' Journal of the Royal Statistical Society Series A, 184(3), 1030-1051,
#' http://doi.org/10.1111/rssa.12699
#' @examples
#' trapdoor
#' trapdoor$plot()
#' @export
trapdoor <- SCM$new("trapdoor",
                    uflist = list(
                      uw = function(n) {return(stats::runif(n))},
                      uz = function(n) {return(stats::runif(n))},
                      ux = function(n) {return(stats::runif(n))},
                      uy = function(n) {return(stats::runif(n))},
                      uwx = function(n) {return(stats::runif(n))},
                      uwy = function(n) {return(stats::runif(n))}
                    ),
                    vflist = list(
                      w = function(uw,uwx,uwy) {
                        return( as.numeric(uw < 0.3 + 0.3*uwx + 0.3*uwy))},
                      z = function(uz,w) {
                        return( as.numeric(uz < 0.4 + 0.4*w))},
                      x = function(ux,z,uwx) {
                        return( as.numeric(ux < 0.3 + 0.5*z - 0.2*uwx ))},
                      y = function(uy,x,uwy) {
                        return( as.numeric(uy < 0.1 + 0.4*x + 0.4*uwy))}
                    )
)

#' SCM "backdoor_md" used in the examples.
#'
#' Variable z fulfills the back-door criterion for P(y|do(x)).
#' Variable z is missing completely at random. The missingness of
#' variables x and y depend on z.
#' @examples
#' backdoor_md
#' backdoor_md$plot()
#' @export
backdoor_md <- SCM$new("backdoor_md",
                       uflist = list(
                         uz = "n : runif(n)",
                         ux = "n : runif(n)",
                         uy = "n : runif(n)",
                         urz = "n : runif(n)",
                         urx = "n : runif(n)",
                         ury = "n : runif(n)"
                       ),
                       vflist = list(
                         z = "uz : as.numeric(uz < 0.4)",
                         x = "ux, z : as.numeric(ux < 0.2 + 0.5*z)",
                         y = "uy, z, x : as.numeric(uy < 0.1 + 0.4*z + 0.4*x)"
                       ),
                       rflist = list(
                         z = "urz : as.numeric( urz < 0.9)",
                         x = "urx, z : as.numeric( (urx + z)/2 < 0.9)",
                         y = "ury, z : as.numeric( (ury + z)/2 < 0.9)"
                       ),
                       rprefix = "r_"
)



#' SCM "credit" used in the credit scoring example.
#'
#' Variable \code{default} is the outcome to be predicted
#' @examples
#' credit$simulate(100)
#' summary(credit$simdata)
#' @export
credit <- SCM$new("credit scoring",
                  uflist = list(
                    u_age = "n : rnorm(n)",
                    u_marital = "n : rnorm(n)",
                    u_gender = "n : rnorm(n)",
                    u_education = "n : rnorm(n)",
                    u_children = "n : rnorm(n)",
                    u_job = "n : rnorm(n)",
                    u_income = "n : rnorm(n)",
                    u_housing = "n : rnorm(n)",
                    u_address = "n : rnorm(n)",
                    u_savings = "n : rnorm(n)",
                    u_credit_history = "n : rnorm(n)",
                    u_length_of_employment = "n : rnorm(n)",
                    u_ethnicity = "n : rnorm(n)",
                    u_credit_amount = "n : rnorm(n)",
                    u_default = "n : rnorm(n)",
                    u1 = "n : rnorm(n)", 
                    u2 = "n : rnorm(n)", 
                    u3 = "n : rnorm(n)", 
                    u4 = "n : rnorm(n)",
                    u5 = "n : rnorm(n)"
                  ),
                  vflist = list(
                    age = "u_age : round(18 + 60*pnorm(u_age)) ",
                    ethnicity = "u_ethnicity : rcateg(c(0.75,0.15,0.10), pnorm(u_ethnicity), softmax = FALSE)",
                    marital = "u_marital, u2, age: 
                                rcateg(cbind(1, 
                                         0.5 + 0.1*(age + 3*u2 - 22), 
                                         0 + 0.2*(age + 2*u2 - 29)), 
                                      pnorm(u_marital), softmax = TRUE)",
                    gender = "u_gender: as.numeric(u_gender > 0)",
                    education = "u_education, u4, u2, gender, age, ethnicity : 
                                  rcateg(cbind(-gender, 
                                  u4 + u2 - gender + 0.1*age + as.numeric(ethnicity == 1),
                                  2*u4 + 2*u2 - 2*gender + 0.2*(age - 22) + 2*as.numeric(ethnicity == 1),
                                  3*u4 + 3*u2 - 3*gender + 0.3*(age - 28) + 3*as.numeric(ethnicity == 1)),
                                   pnorm(u_education), softmax = TRUE, asfactor = FALSE)",
                    children = "u_children, age, ethnicity, education, marital: 
                                 qpois( p = pnorm(u_children),
                                 lambda = pmax(0.2, -0.2 * (ethnicity == 1) + 0.5*(marital == 3) +
                                 as.numeric(age < 45)*(age - 18)/13 + 2*as.numeric(age >= 45)))",
                    job = "u_job, u4, u2, gender, ethnicity, education, children, age:
                            rcateg(cbind( 0.2*(ethnicity != 1) + (education <= 1),
                            u2 + u4 + (education >= 2),
                            u2 + u4 + (education == 3) + 0.01*(age - 28)), 
                            pnorm(u_job), softmax = TRUE)",
                    income = "u_income, u4, job, education, length_of_employment, age:
                              pmax(0, -10000 + 25000*as.numeric(job) + 5000*u4 + 10000*u_income + 2000*education + 
                              200*length_of_employment + 100*(age - 35))",
                    housing = "u_housing, u3, u5, marital, children, age, education, income:
                                rcateg(cbind(7, 0.1*( u3 + u5 + as.numeric(marital) + children + age + education + income/10000)), 
                                pnorm(u_housing), softmax = TRUE)", # 1 rent, 2 own
                    address = "u_address, u5, marital, children, ethnicity, age, income:
                                round( pmin(10,pmax(1, -4 + u_address + u5 + as.numeric(marital) + 
                                      (ethnicity == 1) + age/30 + income/10000 )))",
                    savings = "u_savings, u1, u3, marital, gender, ethnicity, children, education, income, age:
                                pmin(900000 + 10000*exp(u_savings), pmax(0, -5000 + exp(u_savings + u1 + u3 + income/10000 + age/10)))",
                    length_of_employment = "u_length_of_employment, age, education:
                                    pmax(0, pmin(age - 18 - 3*education + u_length_of_employment, 
                                                0.3*(age - 18)  + 3*u_length_of_employment))",
                    credit_amount = "u_credit_amount, u1, income, age, housing, job, savings, marital, children:
                     pmax(5000, 110000 - 4000*abs(age - 40) + 2*(income - 30000) - 20000*as.numeric(housing) + 10000*as.numeric(job) +
                     - 0.2*savings + 20000*as.numeric(marital) + 10000*children + 40000*u_credit_amount + 10000*u1)",
                    default = "u_default, u1, length_of_employment, housing, income, age, savings, education, job, 
                    credit_amount, ethnicity:
                               as.numeric(qlogis(p = pnorm(u_default), 2.2*income + 10000*education + 
                               10000*as.numeric(job) + 3000*length_of_employment + 10000*u1 - 0.7*credit_amount +
                                5000*as.numeric(ethnicity), scale = 5000) < 0)" # 1 defaults
                  )
)
