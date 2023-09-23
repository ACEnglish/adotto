
paths.txt holds metadata about each of the benchmarking results

We want to look at these results and use them to validate the benchmark.

1) How much of a difference does refine make on the results?
2) How informative are the reports?
    - Can we use variant / region / laytr to get a better understanding of what any particular tool is doing?
    - Want to bring attention to the subsets stratifications.
3) Are there sites inaccessible from one technology / technique vs another?
    - I bet there are sites where non-nist WGS all miss but every TR caller captures
    - I bet there are sites long reads capture but none of the short reads capture
4) How many locations are universally FP or FN?
    - These may be indicative of low quality sites to remove from the benchmark.

Okay, I've made the data frames in each of the results.
And I have all.summary.jl with the different stats for truvari.
And I made all the laytr reports
