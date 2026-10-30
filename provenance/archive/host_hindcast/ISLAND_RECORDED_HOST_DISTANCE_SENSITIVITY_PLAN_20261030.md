# Exploratory source-distance tilt sensitivity (2026-10-30)

**Follows:** [the frozen uniform recorded-host comparison](ISLAND_RECORDED_HOST_EXACT_PLAN_20261030.md), whose 20 observed recorded-host cells, exact uniform p=0.04410064, and three-block leave-one-out instability have already been examined. This is a separately labelled *post-result model-dependent sensitivity*, **not** a replacement for the locked uniform finite-law result and **not** preregistration.

## Question

How much of the weak documented early-host co-recording signal may arise from the pronounced proximity of later butterfly occurrence cells to previously recorded butterfly source islands?

## Locked alternative conditional law and calculations

Use exactly the same 2,385 risk cells, 66 observed positives, 12 block allowed masks, fixed species×longitude×latitude-zone row totals, and island column totals. All original host-taxonomy statuses are unchanged. For each complete permissible binary table T, let D(T) be the sum of log(1 + nearest previously *recorded-island* coast gap in km) over later-positive cells, as in the source-distance paper. Let H(T) be the number of documented prior catalogue-compatible host recording cells among those positives.

Define a one-parameter **illustrative source-distance tilt** on this exact finite set:

P_lambda(T) proportional to exp(−lambda × D(T)), lambda >= 0.

Every retained table receives its normalized weight; enumerate each block completely, with no Monte Carlo, and convolve its weighted integer-H distribution across blocks. The distance-tilt parameter lambda* is fixed by the moment-matching equation E_lambda[D] = D_observed (when a finite nonnegative root exists). Solve by bracketed bisection to tight numerical tolerance, without using host scores to fit lambda. Also show lambda=0 (the earlier uniform null), lambda*/2, and 2×lambda* as illustration; do not selectively choose whichever gives desired host outcome.

Primary diagnostic: normalized weighted P_lambda*(H >= 20) and E_lambda*[H]. **This is a plug-in descriptive tail conditional on fitted distance preference, not an exact p-value for a composite nuisance-parameter null.** This model is also not a demonstrated real dispersal probability; closest indexed source need not be a source population, and recording effort remains.

## Checks

- At lambda=0, reobtain uniform H expectation 17.06018423356553 and inclusive tail 6,487,980 / 147,117,600.
- Each block's full table list reproduces existing exact cardinalities 52,542, 10, 28, 10 or 1.
- Total D observed and lambda=0 expected D reproduce the established independent enumeration.
- Total normalized weighted H probabilities sum to 1 for each lambda, and the fitted E_lambda*[D] matches observed D to numerical tolerance.
- Any change of the ecological claim remains prohibited. Even a large source-distance-adjusted descriptive tail would not independently identify host quality or causal host availability.
