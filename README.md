
# SIR-Hawkes: Linking Epidemic Models and Hawkes Processes to Model Diffusions in Finite Populations
This repository contains:
 - Scripts for SIR-Hawkes project: simulation, modeling with SIR and HawkesN.
 - Three datasets consist of tweet cascades.
 - A hands-on tutorial to walk you through some main components of the project: simulation, modeling and population size distribution.

### Citation
The project was introduced in this [paper](https://arxiv.org/abs/1711.01679):
```
Marian-Andrei Rizoiu, Swapnil Mishra, Quyu Kong, Mark Carman, Lexing Xie. 2018. SIR-Hawkes: 
Linking Epidemic Models and Hawkes Processes to Model Diffusions in Finite Populations. . 
In WWW 2018: The 2018 Web Conference, April 23–27, 2018, Lyon, France. ACM, New York, NY, 
USA, 10 pages. https://doi.org/10.1145/3178876.3186108
```
### License
Both dataset and code are distributed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license, a copy of which can be obtained following this link. If you require a different license, please contact us at Marian-Andrei@rizoiu.eu or Lexing.Xie@anu.edu.au.

# SIR-Hawkes tutorial

### required packages:
    - nloptr
    - parallel
    - data.table

### 1. Preliminary
We need to first load all required packages for simulation and modeling cascades.


```R
library(parallel)
source('scripts/functions-SIR-HawkesN.R')
source('scripts/functions-size-distribution.R')
```

### 2. Stochachastic R simulation
We then simulate 20 stochastic SIR realizations. In this step, we chose a set of parameters (<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/e84961284909358d9ae8cb817dd85569.svg?invert_in_darkmode" align=middle width=245.536005pt height=22.83138pt/>) for simulation. Given those simulated events, we are going to fit them with both `SIR` model and our proposed `HawkesN` model to see their modeling performance.


```R
params.S <- c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1)
nsim <- 20
simdat <- replicate(
    n = nsim,
    generate.stochastic.sir(params = params.S, Tmax = 11, hide.output = T)    
)
```

Let's take a look at the simulated data (only the first 20 events of the first simulation were shown).

One simulation is identified as four components:
    - relative times
    - susceptible population size at each time
    - infected population size at each time
    - recovered population size at each time


```R
as.data.frame(simdat[,1])[1:20,]
```


<table>
<thead><tr><th scope=col>time</th><th scope=col>S</th><th scope=col>I</th><th scope=col>R</th><th scope=col>C</th></tr></thead>
<tbody>
	<tr><td>0.000000000</td><td>1000       </td><td>300        </td><td>0          </td><td>300        </td></tr>
	<tr><td>0.001861149</td><td> 999       </td><td>301        </td><td>0          </td><td>301        </td></tr>
	<tr><td>0.002054043</td><td> 998       </td><td>302        </td><td>0          </td><td>302        </td></tr>
	<tr><td>0.002297054</td><td> 997       </td><td>303        </td><td>0          </td><td>303        </td></tr>
	<tr><td>0.003679539</td><td> 996       </td><td>304        </td><td>0          </td><td>304        </td></tr>
	<tr><td>0.005308822</td><td> 995       </td><td>305        </td><td>0          </td><td>305        </td></tr>
	<tr><td>0.010011856</td><td> 995       </td><td>304        </td><td>1          </td><td>305        </td></tr>
	<tr><td>0.013833165</td><td> 995       </td><td>303        </td><td>2          </td><td>305        </td></tr>
	<tr><td>0.015425608</td><td> 994       </td><td>304        </td><td>2          </td><td>306        </td></tr>
	<tr><td>0.018386120</td><td> 993       </td><td>305        </td><td>2          </td><td>307        </td></tr>
	<tr><td>0.018890458</td><td> 992       </td><td>306        </td><td>2          </td><td>308        </td></tr>
	<tr><td>0.019177649</td><td> 992       </td><td>305        </td><td>3          </td><td>308        </td></tr>
	<tr><td>0.032191757</td><td> 991       </td><td>306        </td><td>3          </td><td>309        </td></tr>
	<tr><td>0.035561494</td><td> 990       </td><td>307        </td><td>3          </td><td>310        </td></tr>
	<tr><td>0.037340979</td><td> 989       </td><td>308        </td><td>3          </td><td>311        </td></tr>
	<tr><td>0.038834789</td><td> 988       </td><td>309        </td><td>3          </td><td>312        </td></tr>
	<tr><td>0.041421224</td><td> 987       </td><td>310        </td><td>3          </td><td>313        </td></tr>
	<tr><td>0.041912680</td><td> 986       </td><td>311        </td><td>3          </td><td>314        </td></tr>
	<tr><td>0.043324751</td><td> 985       </td><td>312        </td><td>3          </td><td>315        </td></tr>
	<tr><td>0.051393662</td><td> 984       </td><td>313        </td><td>3          </td><td>316        </td></tr>
</tbody>
</table>



### 3. Fit Stochastic SIR on simulated cascades
We fit stochastic SIR in the following steps:
 - Choose a starting point for all parameters.
 - Apply LBFGS algorithm for optimizing the likelihood function of `SIR` model (This step might take quite a lot of time).


```R
# initial fitting point for each execution
params.fit.start <- c(N = 0.1, I.0 = 0.1, gamma = 0.1, beta = 0.1)

.cl <- makeCluster(spec = min(nsim, detectCores()), type = 'FORK')
results <- parSapply(cl = .cl, X = 1:nsim, FUN = function(i) {
    mysim <- as.data.frame(simdat[, i])
    return(fit.stochastic.sir(mysim, params.fit.start))
})
stopCluster(.cl)

# reconstruct result data format
res <- as.data.frame(results[1,])
names(res) <- 1:nsim             
res <- as.data.frame(t(res))     
res<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/9914b698d46caf7a7b1b6120c64347ad.svg?invert_in_darkmode" align=middle width=1244.40855pt height=1268.8599pt/>\{\beta, \gamma, N\}<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/5d211b0e3e6971e44cb61ea68560507b.svg?invert_in_darkmode" align=middle width=233.832555pt height=22.83138pt/>\lambda^I(t)<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/2c2b722c19a50b7108d6c84c29eb9c70.svg?invert_in_darkmode" align=middle width=436.804005pt height=22.83138pt/>\{\mu, \kappa, \theta, N\}<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/734ae0cd7c1abb8d83c3326e1d3c141c.svg?invert_in_darkmode" align=middle width=124.374855pt height=22.83138pt/>\lambda^H(t)<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/503de8ef1a6ab1b251dd7a66a6ddcad1.svg?invert_in_darkmode" align=middle width=29.343765pt height=22.46574pt/>\mathcal{T} = \{\tau_1, \tau_2, \ldots\}<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/53b23f78d233638ddd9daab1911ceaee.svg?invert_in_darkmode" align=middle width=588.415905pt height=22.83138pt/>\lambda^I(t)<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/672817601eea255514a955612d495cd9.svg?invert_in_darkmode" align=middle width=32.053065pt height=14.15535pt/>\mathcal{T}<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/3b874d4cad3fc06b8995b77a4ef61ee1.svg?invert_in_darkmode" align=middle width=52.278765pt height=22.83138pt/>\lambda^H(t)<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/f11ece218135bb1a349b13505cb763d3.svg?invert_in_darkmode" align=middle width=182.261805pt height=27.65697pt/>\mu = 0<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662925pt height=14.15535pt/>\beta = \kappa \theta<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/24ee684c2922b0d32c54a34089c92ec0.svg?invert_in_darkmode" align=middle width=4.5662925pt height=14.15535pt/>\gamma = \theta<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/6e11b99d0336defa4a4545df7dc453ea.svg?invert_in_darkmode" align=middle width=700.27485pt height=118.35615pt/>gamma <- res<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/5697e373afedc5962ce6f6b59694a896.svg?invert_in_darkmode" align=middle width=60.91932pt height=22.83138pt/>beta <- res<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/8b302efaf0737751585d22c22c415bbf.svg?invert_in_darkmode" align=middle width=53.89461pt height=22.46574pt/>theta
prnt <- rbind(params.S[c('N', 'gamma', 'beta')], 
              apply(X = res[, c('N', 'gamma', 'beta')], MARGIN = 2, FUN = mean, na.rm = T),
              apply(X = res[, c('N', 'gamma', 'beta')], MARGIN = 2, FUN = sd, na.rm = T))
rownames(prnt) <- c('theoretical', 'median', 'sd')
print(prnt[, c('N', 'gamma', 'beta')], digits = 2)
```

                   N gamma beta
    theoretical 1300   0.2  1.0
    median      1295   3.5  6.4
    sd            11   8.5 13.7


### 5. Plotting size distribution
In this section, we will study the probability distribution of population size for a given set of parameters. 


```R
# theoretical parameters shown previously
params.S <- c(N = 1300, I.0 = 300, gamma = 0.2, beta = 1)
params.H <- c(K = 5, c = 0.001, theta = 0.2, N = 1300)
```

Let's calculate the transition matrix and the states mapping matrix for an SIR model with the given parameters.


```R
.transition <- construct.transition.matrix.and.states(params = params.S)
```

Then we construct the event history in required format for HawkesN and we choose the first cascade from our simulated data.


```R
seen_perc <- 0.5
history <- as.data.frame(simhistory[,1])
seenEvents <- round(nrow(history) * seen_perc)
```

Then we get two size distributions for two scenarios:
 - when only 1 event is observed (apriori probability size distribution)
 - when half of events are observed (aposteriori probability size distribution)


```R
# compute size probs at event 1 and current event
size.est.at.zero <- get.size.distribution(params = params.H, .transition = .transition)
size.est.at.end <- get.size.distribution(params = params.H, .transition = .transition, history = history[seq(seenEvents),])
```

We then plot both apriori probability size distribution and aposteriori probability size distribution.


```R
# plot our cascade
matplot(cbind(size.est.at.zero<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/6be94403d6d6f5823e7d0115e22125cc.svg?invert_in_darkmode" align=middle width=192.145305pt height=22.83138pt/>final.state), 
        col = c('black', 'blue'), type = 'l', log = 'y', lty = c(1, 1), lwd = 3,
        xlab = 'Cascade final size', ylab = 'Probability',
        main = sprintf('Probability distribution of cascade final size\nfitted on %d seen events (N = %.2f)\n(seen %%: %.2f)', 
                       seenEvents, params.H['N'], seen_perc) )
abline(v = seenEvents, lty = 3, col = 'gray40')
abline(v = size.est.at.end<img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/9bc9462d58d1b876835015ec4910f794.svg?invert_in_darkmode" align=middle width=926.04435pt height=47.67114pt/>theo.mean), 
                   sprintf('Observed final size (%d)', nrow(history)) ), 
        lty = c(1, 1, 3, 2, 1), lwd = c(3, 3, 1, 1, 1), col = c('black', 'blue', 'gray40', 'darkmagenta', 'red'), bty = 'n')
```

    Warning message in xy.coords(x, y, xlabel, ylabel, log = log):
    “642 y values <= 0 omitted from logarithmic plot”Warning message in xy.coords(x, y, xlabel, ylabel, log):
    “1 y value <= 0 omitted from logarithmic plot”


![png](images/output_29_1.png)


Several observations (Sec 6.3 in our paper):
 - The apriori probability size distribution shows two maxima. This provides the following explanation for the general perceived unpredictability of online popularity. For cascades showing a bi-modal apriori size distribution, there are two likely outcomes: either it dies out early or it reaches a large size compared to the maximum population <img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/f9c4988898e7f532b9f826a75014ed3c.svg?invert_in_darkmode" align=middle width=14.999985pt height=22.46574pt/>. At time <img src="https://rawgit.com/in	git@github.com:computationalmedia/sir-hawkes/master/svgs/477a717e18587a5e8605780ca167c322.svg?invert_in_darkmode" align=middle width=36.07296pt height=21.18732pt/> is it impossible to di erentiate between the two outcomes.
 - The aposteriori probability distribution reflects the information gained from the observed events and it shows a single maximum towards the higher size values. The more events we observe, the higher the likelihood of the true value of cascade size.
