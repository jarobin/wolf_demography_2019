# Script to make and submit SLIM job for wolf simulations in Robinson et al. 2019.
# Model is based on population sizes and divergence times from Fan et al. 2016,
# DOI:10.1101/gr.197517.115

# USAGE: ./make_submit_slim_wolf.job.sh [u] [h] [t]
# Example: ./make_submit_slim_wolf.job.sh 1e-8 0.5 450000


# Specify path to SLIM executable
SLIM=/opt/SLiM/bin/slim

# Set u, the per-site mutation rate
u=${1}

# Set h, the dominance coefficient for deleterious mutations 
# (e.g. h=0.0 for fully recessive, h=0.5 for additive)
h=${2}

# Set t, the number of burn-in generations to build up variation in the ancestral 
# population. Recommended value is 10N, where N is the number of individuals in the 
# ancestral population.
t=${3}


# Make script
cat > slim_wolf_1000genes_u${u}_h${h}_${t}burn.job << EOM

// Wolf Sim
// slim_wolf_1000genes_u${u}_h${h}_${t}burn.job

// INITIALIZE

initialize() {
initializeMutationRate(${u}); 

// Distribution of selection coefficients for nonsynonymous mutations from 
// Kim et al. 2017, DOI:10.1534/genetics.116.197145
initializeMutationType("m1", ${h}, "g", -0.01314833, 0.186); 

// Neutral mutations
initializeMutationType("m2", 0.5, "f", 0.0);

// Synonymous:nonsynonymous mutation ratio of 1:2.31 from Huber et al. 2017,
// DOI:10.1073/pnas.1619508114
initializeGenomicElementType("g1", c(m1,m2), c(2.31,1.0));

for (i in 1:1000){
	initializeGenomicElement(g1, ((i-1)*1000)+(i-1), (i*1000)+(i-2) );
}

// Number of genes per chromosome for 1000 genes, assuming genes are evenly distributed
// and chromosome lengths correspond to relative lengths in the dog genome.
gene_nums=c(56,39,42,40,40,35,37,34,28,31,34,33,29,28,29,27,29,25,24,26,23,28,24,21,23,18,21,19,19,18,18,18,14,19,12,14,14,11);

// Set recombination rates:
// Recombination rate in between "chrosomomes" is 50% (i.e. unlinked).
// Recombination rate within genes is 0.
// Recombination rate in between genes is 1e-3, the effective recombination rate for 
// 100 kb noncoding sequence with a per-site recombination rate of 1e-8 in between 
// each gene.
rates=NULL;
for (i in 1:(size(gene_nums)-1)){
	rates=c(rates, 0, rep(c(1e-3, 0), gene_nums[i-1]-1), 0.5);
}
rates=c(rates, 0, rep(c(1e-3, 0), gene_nums[size(gene_nums)-1]-1));

ends=NULL;
for (i in 1:${g}){
	ends=c(ends, (i*1000)+(i-2), (i*1000)+(i-1));
}
ends=ends[0:(size(ends)-2)];

initializeRecombinationRate(rates, ends);
}


// SIMULATE

// Begin the simulation with burn-in period for ancestral population.
1 { 
	sim.addSubpop("p1", 45000); 
} 

// After burn-in, establish two populations: one large (N. American wolves) and one small
// (Tibetan wolves).
${t} {
	sim.addSubpopSplit("p2", 2500, p1);
	p1.setSubpopulationSize(17350);
}

// Keep track of generation number in log file every 1,000 generations 
// (overwrites itself).
1:$((${t} + 4167)) late() {
	if (sim.generation % 1000 == 0){
		writeFile("slim_wolf_1000genes_u${u}_h${h}_${t}burn.gen", paste(sim.generation));
	}
}

// After 4,167 generations the simulation ends and output is written to stdout
$((${t} + 4167)) late() {
  // Output population sizes for easy reference
  cat("#p1: " + size(p1.individuals) + " individuals; p2: " + size(p2.individuals) + " individuals.\n");
  
  // File header:
  // Mutation id,
  // Mutation type,
  // Selection coefficient,
  // Age of mutation in generations,
  // Subpopulation it arose in,
  // Number of heterozygote derived in p1,
  // Number of homozygote derived in p1,
  // Number of heterozygote derived in p2,
  // Number of homozygote derived in p2
  // Note: these are genotype counts, not allele counts
  cat("mutid,type,s,age,subpop,p1numhet,p1numhom,p2numhet,p2numhom\n");
  
  // Collect stats for each mutation in the simulation
  for (mut in sim.mutations){
    id = mut.id;
    s = mut.selectionCoeff;
    subpop = mut.subpopID;
    age = sim.generation - mut.originGeneration;
    type = mut.mutationType;
    if (type == m1){
      type = "m1";
    } else if (type == m2){
      type = "m2";
    }

    // Initialize genotype counts
    p1numhet = 0;
    p1numhom = 0;
    p2numhet = 0;
    p2numhom = 0;
    
    // Count number of derived heterozygotes and homozygotes in p1
    for (p1i in p1.individuals){
      gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
      if (gt == 1){
        p1numhet = p1numhet + 1;
      } else if (gt == 2){
        p1numhom = p1numhom + 1;
      }
    }
    
    // Count number of derived heterozygotes and homozygotes in p2
    for (p2i in p2.individuals){
      gt = sum(c(p2i.genomes[1].containsMutations(mut), p2i.genomes[0].containsMutations(mut)));
      if (gt == 1){
        p2numhet = p2numhet + 1;
      } else if (gt == 2){
        p2numhom = p2numhom + 1;
      }
    }

    // Print results
    cat(id + ",");
    cat(type);
    cat("," + s + "," + age + "," + subpop + "," + p1numhet + "," + p1numhom + "," + p2numhet + "," + p2numhom + "\n");
  }
}

EOM


# Execute script

${SLIM} ./slim_wolf_1000genes_u${u}_h${h}_${t}burn.job
