string helpstring = R"(The following are possible flags for assemblyCpp:

-runTime=x : sets runTime to x miliseconds in windows and x microseconds in linux, default is unlimited

-enumMax=x: sets the number of subgraphs the algorithm is allowed to enumerate to x, default is 50,000,000. Setting this too large will cause memory issues for big graphs

-pathwayFolder=x: sets the pathway folder for the string assembly algorithm to x

-pathway=x: if x is 1 (the default), a pathway file is generated, else if x is 0 no pathway will be generated

-removeHydrogens=x: AssemblyCpp removes explicit hydrogens by default. If x is set to 0, the calculation takes into account explicit hydrogens. Only applies to molfile inputs.

-disjointCompensation=x: Compensates for joint assembly spaces if x=1 by subtracting 1 for each disjoint graph after the first, else does nothing if x=0 (default)

-memTest=x: Turn on memory usage tracking if x=1 (linux only), else does nothing if x=0 or windows (default)

-writeIntermediateMAs=x: Write out the smallest MA found so far and the time it was found at into a separate file if x=1, else does nothing if x=0 (default)
)";

/**
 * @brief Produces helpstring
 * 
 */
void help()
{
    cout << helpstring << endl;
}