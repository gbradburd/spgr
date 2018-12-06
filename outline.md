# Things to Measure

1. Where are the things?

    a. *Population density:* 

        - How many individuals live in a given area? in total?

        - How has that number of individuals fluctuated over time? (fast changes as opposed to slow changes)

    b. *Fecundity/reproduction:*

        - How many new offspring are produced per year in a given area? in total?

        - How do areas contribute in the *long term*?

        - picture: 
            in an example with higher fecundity in some areas 
            (gradient of mean fecundity superimposed on a bumpy landscape?)
            and small bias in dispersal towards lower fecundity area

            * map of population density juxtaposed with
            * a map of per capita number of offspring produced
            * and a map of long-term fitness 

    c. Related statistics:

        - on-the-ground census
        - $N_e$


2. How do they move?

    a. *Dispersal* (individual, diffusive movement; $\sigma$):

        - What is the mean distance between parent and child birth location (in a homogeneous area)?
        - Is there a mean bias in parent-child displacement?

    b. *Population flux:*

        - How many individuals are born each year whose parents were born on the other side of a given boundary?
            (separately for each direction)

            * picture: adjacent valleys with lower pop density along the intervening ridge,
                with arrows denoting all parent-child relationships that cross the ridge 

            * picture: added to the above, a line drawn down the middle of a valley, with arrows crossing as before

            * relationship to dispersal: birth rate times population density times length of boundary times $\sigma/\sqrt{2\pi}$ in each direction 
                if dispersal is Gaussian and boundary is straight

    c. Related statistics:

        - $m$ between "discrete pops"
        - resistance distance


3. How different are they in different places?

    a. *(Relative) genetic differentiation:*

        - How much is genetic diversity partitioned across geography?
            How does the expected genetic divergence (or, TMRCA) for nearby individuals compare to distant ones?
            
            * $F_{ST}(x,y) = 1 - (\pi(x,x) + \pi(y,y))/(2 (\pi(x,x) + 2 \pi(x,y) + \pi(y,y)))$, 
                i.e., defining mean $F_{ST}$ at distance $x$ to be one minus ratio of mean TMRCA at each location
                to mean TMRCA between individuals from the pops pooled

            * relate to the classic IBD curve

        - How much local inbreeding is there?

            * $\pi(x,x)$? is not heterozygosity, $\pi(x,x)$ should be $\lim_{y \to x} \pi(x,y)$,
                so we can measure this as heterozygosity minus $\pi(x,x)$

            * picture: map of heterozygosity across space in a habitat with one big valley surrounded by some small valleys


    b. *Geographic distribution of ancestry* (and the opportunity for local adaptation):

        - How much of the average individual's ancestry is still within a given region, $T$ years in the past?

        - How far back in time do you have to go before the average proportion of ancestry that an individual inherits
            from a given region drops below (say) 90%?
            
            * picture: individuals shaded according to the proportion of ancestry of a given small cluster of samples
                they contribute, at several points back in time


    c. Related statistics:

        - $F_{ST}$, but already covered this

        - relative positions on a PC plot?


4. How have things changed over time?

    a. Scenario: constant range but more-or-less gradually changing population density due to changing balance of competition/resource availability.

        - Relate to $N_e$-over-time things.

    b. Scenario: postglacial expansion and secondary contact.

        - How much does each individual inherit from each glacial refugium?

            * picture: pie charts of this quantity.

        - Relate to admixture.

    c. Scenario: orogeny induces phylogeny.

        - What proportion of ancestry of one side of the mountain comes from the other side of the mountain,
            as a function of time ago?

            * picture: across a barrier with gradually increasing strength,
                plot of this quantity, compared to if the mountain had not grown.

    d. Scenario: range expansion and exponential growth.

        - Where did the range expansion originate?

        - How fast did it spread?

        - How much surfing was there on the wave front?

