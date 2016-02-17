
require(phylosim)


# Program = the program for which the output is formatted for
# tree = tree
# length = number of characters to be simulated
# replicates = number of replicates to be simulated
# file = the naming convention required (eg. 100_char_Rate_1 etc.)
# directory = the directory that the files are stored in for the HPC, leave as default if program=PAUP
# pi = a vector of length 4, which sums to 1, this is the equilibrium freqs (T,C,A,G) respectively.
# kappa = the number of transtions relative to transversions
# rate.distro= the distribution used to measure site specific rate heterogeneity, string of 'Gamma' or 'Lognormal'
# parameter = the variance of the lognormal distribution with arthimetric mean = 1, or the shape of the gamma distribution
sim.and.write <- function(tree, length, replicates, file="OUTPUT", directory="/foo/bar/cluster", program='Both', pi=c(0.25,0.25,0.25,0.25), kappa=1, rate.distro='Gamma', parameter=1, samples){

	PSIM_FAST <- TRUE # This enables the PhyloSim 'Fast and Messy' Mode

	original <- getwd() # Get the current dir
	directory.creation(filename=file, replicates, samples) #Make a home for each replicate

	setwd(original) # go home

	p <- HKY(rate.params=list("Alpha"=kappa, "Beta"=1), base.freqs=pi) # Form the HKY Model

	# 							Some HKY Notes                                                                 #
	########################################################################################################################################
	# Alpha = Transitions
	# Beta = Transversions
	# Kappa = Transitions/Transversions
	# Common practice is to set Beta to 1 and to allow Alpha to vary
	# This still allows kappa to reperesent the ratio
	# Transitions occur at a Higher frequency than Transversions, ie. Alpha>Beta
	# The rate matrix therefore appears as: http://upload.wikimedia.org/math/0/4/b/04b4d3860f67108ce9105d218def97b5.png
	# Where Kappa represents a Transition, and Beta=1=Transversion rate
	# if Kappa=1 and Pi = [0.25,...,0.25] we collapse to JC69 (p<-HKY(rate.params=list( "Alpha"=1,"Beta"=1),base.freqs=c(0.25,0.25,0.25,0.25));plot(p))
	# Importantly, Under this parameterisation we change a single value, Alpha (and therefore Kappa) instead of 6+ in GTR
	########################################################################################################################################
    
	root.seq <- NucleotideSequence(length=length, processes=list( list(p) ))# Get an ancestral sequence of length=length made up exclusively of '?'s

	underlying_rates <- runif(0.1, 10, n = samples) # Take n samples from a uniform distribution

	for(m in 1:samples){

		setwd(file.path(original,paste0(file, '/', form.filename(file, m), sep = ""))) # Get into the correct dir for the particular sample
		sampledir <- getwd() # store the sample dir in a variable
	
		write.table(underlying_rates, "Rate_Multipliers_For_BrLen_Alteration.txt") # Write the n samples from the uniform distribution to file	
		print(paste("Rate", m, "of", samples, "=", underlying_rates[m]))


		tree_m <- tree # Get a tree for each sample
		tree_m$edge.length <- tree_m$edge.length*underlying_rates[m] # set its brlens to be original*sample[m]

		for (j in 1:replicates){ # For each replicate do the simulations

			clearStates(root.seq) # Clear any ancestral states

			if (rate.distro == 'Gamma'){
				plusGamma(root.seq, p, parameter, ncat = "cont") # Sample a rate for each site from a continuous Gamma
				getRateMultipliers(root.seq,p);
			}

			if (rate.distro == 'Lognormal')	{ # The mean in real space is set to 1, as it is for Gamma when a=b, we therefore use the arithemtic moments equations to set up our parameters such that arth mean=1 and arth variance is a free parameter
				      

	# 							Some Lognormal Distribution Notes                                              #
	########################################################################################################################################
	# The Lnorm distro is paramtised as mean and sd, in log space
	# We need to set values in arithemtic, real space
	# http://upload.wikimedia.org/math/0/e/d/0ed02e0e945ce2bf85d23f0869a034f4.png
	# These Equations will allow us to set mean =1 whilst varying the variance in arithemtic space
	# Setting up equations to automate the conversion of these values into meanlog and sdlog seems senisble
	########################################################################################################################################


				mu <- log(1)-(0.5*(log(1+(parameter/1^2)))) # Calculate mu for arithmetic mean=1 and specified var
				sigma <- log(1+(parameter/1^2)) #calculate sig2 for arMean=1 and specified var
				sigma <- sqrt(sigma) # get the root, leaving sigma, the parameter required

				multipliers <- rlnorm(n=length,meanlog=mu,sdlog=sigma) # Use the rlnorm function to take n=length random draws from a lognormal distribution
				setParameterAtSites(root.seq, p, "rate.multiplier", multipliers) #Use the setParamterAtSites function to override the need for a gamma distribution
			}	
		sampleStates(root.seq)	# Get a new ancestral sequence for each replicate by using the state freqs

		print(j) # Which replicate is it?

		setwd(file.path(sampledir, paste0(file, "_", m, "_", j))) # Get into the particular replicate dir for the particular sample 
		#setwd(file.path(original,paste0(file,'/',form.filename(file, j), sep = ""))) # Get into the correct dir

		output <- seq.sim(tree=tree_m, length=length, root.seq)# Do the simulation

		if (program == 'MrBayes'){
			write.nexus.MrBayes(filename=paste(filename = form.filename(file,j),".NEX",sep = ""), table=as.matrix(output)) # Write it to file
			setwd(file.path(original, file, "Qsubs")) # Go to the Qsub folder
			form.shellscript(filename = file, scriptname = paste(form.filename(file, j),".sh",sep = ""),directory=directory,iteration=j) # Write a PBS script
			setwd(file.path(original, file)) # Go Home
			} 

		if (program == 'PAUP'){ 	# We will ignore the directory argument
			write.nexus.PAUP(filename = paste(filename = form.filename(file, j), ".NEX", sep = ""), table=as.matrix(output), tree= tree) # Write it to file			
			setwd(file.path(original, file)) # Go Home
			}

		if (program =='Both'){ 	# We will write two output files and a qsub for MrBayes
			write.nexus.PAUP(filename = paste(filename = form.filename(file, m), "_", j, "_PAUP.NEX", sep = ""), table=as.matrix(output), tree = tree) # Write it to file
			write.nexus.PAUP.nouninformative(filename = paste(filename = form.filename(file, m), "_", j, "_PAUP_NoInvAut.NEX", sep = ""), table = as.matrix(output), tree = tree) # Write it to file
			write.nexus.MrBayes(filename = paste(filename = form.filename(file, m), "_", j, "_MB.NEX", sep = ""), table=as.matrix(output)) # Write it to file
			setwd(file.path(original, file, "Qsubs")) # Go to the Qsub folder
			form.shellscript(filename = file, scriptname = paste(form.filename(file, m), "_", j, ".sh", sep = ""),directory = directory, iteration = j, sample = m) # Write a PBS script
			setwd(file.path(original, file)) # Go Home
			}
	}
}
	if (program == 'PAUP'){
			unlink("Qsubs", recursive=TRUE) # If PAUP we dont need the shell script dir
	}

	setwd(original) # go home	
}


seq.sim <- function(length, tree, root.seq){
    
    sim <- PhyloSim(root.seq = root.seq, phylo = tree) #Create the simulaton object
    
    Simulate(sim) #Do the simultion
    
    z <- sim$alignment
    
    nodes<-setdiff(row.names(z), tree$tip.label) #get which rownames are internal

	node.row <- c()
	for(i in 1:length(nodes)){ #Now get the actual rows from the names
   		 node.row[i] <- which(row.names(z) == nodes[i])
	}
    z <- z[-node.row,] # remove the rows containing internal values

    z <- gsub("T", "0", z) #recode for Purine/Pyrimidine
    z <- gsub("A", "1", z)
    z <- gsub("C", "0", z)
    z <- gsub("G", "1", z)
 
    colnames(z) <- NULL # Remove the colnames
    
 return(z)
}


# This Function writes the simulated table to a NEXUS formatted file
# For use with MrBayes
# Filename is the name of the output file to be created
# Table is the alignment created by PhyloSim
write.nexus.MrBayes <- function(filename, table){
    

    #Open a connection and write the file
    #This is based on the write.nexus.data.R function in "APE"
    zz<-file(filename,"w")
    "fcat" <- function (..., file = zz) #Function to make concatenation easier
    {
        cat(..., file = filename, sep = "", append = TRUE)
    }
    indent<- "  "
    ntaxa<-length(table[,1])
    nchar<-length(table[1,])
    fcat("#NEXUS\n\n\n")
    NCHAR<- paste("NCHAR=", nchar, sep="")
    NTAX <- paste("NTAX=", ntaxa, sep = "")
    fcat("BEGIN DATA;\n")
    fcat("DIMENSIONS  ",NTAX,indent,NCHAR,";\n")
    fcat("FORMAT DATATYPE=STANDARD;\n")
    fcat("MATRIX\n\n\n\n\n\n")

    
    for (i in 1:length(table[,1])){    	
        fcat("\n", row.names(table)[i], indent, table[i,], "\n\n") 
    }
    
    
    fcat(indent, ";\n\n\n\n\n")
    fcat("begin mrbayes;\n\n\n")
    fcat("set autoclose=yes nowarn=yes;\n")
    fcat("lset coding=all;\n\n\n")
    fcat("lset rates=gamma;\n\n\n") # using Mk not Mkv
    fcat("lset Ngammacat=4;\n\n\n") #4 gamma cats, the default
 
    
    fcat("[mcmc settings]\n\n")
    fcat("mcmc temp=0.1 nchain=4 samplefreq=50 printfr=100 nruns=2 ngen=500000 mcmcdiagn=YES stoprule=yes stopval=0.01;\n\n\n") #This gives us 10,000 samples and stops runs when they get to
    fcat("sumt nruns=2;\n")
    fcat("sump nruns=2;\n")
    fcat("END;\n")
    close(zz) #Close the connection
    

}

# This Function writes the simulated table to a NEXUS formatted file
# For use with PAUP
# Filename is the name of the output file to be created
# Table is the alignment created by PhyloSim
# tree is the tree for which simulations are conducted on, we will write it in the PAUP block
write.nexus.PAUP <- function(filename, table, tree){
    
    tree.to.write<-write.tree(tree) #get the tree as a newick formatted object
    
    
    #Open a connection and write the file
    #This is based on the write.nexus.data.R function in "APE"
    zz<-file(filename,"w")
    "fcat" <- function (..., file = zz) #Function to make concatenation easier
    {
        cat(..., file = filename, sep = "", append = TRUE)
    }
    indent<- "  "
    ntaxa<-length(table[,1])
    nchar<-length(table[1,])
    fcat("#NEXUS\n\n\n")
    NCHAR<- paste("NCHAR=", nchar, sep="")
    NTAX <- paste("NTAX=", ntaxa, sep = "")
    fcat("BEGIN DATA;\n")
    fcat("DIMENSIONS  ",NTAX,indent,NCHAR,";\n")
    fcat("FORMAT DATATYPE=STANDARD;\n")
    fcat("MATRIX\n\n\n\n\n\n")

    
    for (i in 1:length(table[,1])){    	
        fcat("\n", row.names(table)[i], indent, table[i,], "\n\n") 
    }
    
    
    fcat(indent, ";\n\n\n\n\n")
    fcat("END;\n")

    fcat("begin trees;\n\n\n")
    fcat("tree tree_1 = [&U]  ",tree.to.write,"\n\n")
    fcat(";\n\n")
    fcat("END;\n")


    fcat("begin paup;\n\n\n")

	fcat("set increase=auto autoinc=100;\n")
	fcat("hs nreps=1000 swap=tbr addseq=random;\n")
	fcat("savetrees file=equal.tre;\n\n")

	fcat("PSet goloboff=yes gk=2;\n") 
	fcat("hs nreps=1000 swap=tbr addseq=random;\n")
	fcat("savetrees file=gk2.tre;\n\n")

	fcat("PSet goloboff=yes gk=3;\n") 
	fcat("hs nreps=1000 swap=tbr addseq=random;\n")
	fcat("savetrees file=gk3.tre;\n\n")

	fcat("PSet goloboff=yes gk=6;\n") 
	fcat("hs nreps=1000 swap=tbr addseq=random;\n")
	fcat("savetrees file=gk6.tre;\n\n")

	fcat("PSet goloboff=yes gk=10;\n")
	fcat("hs nreps=1000 swap=tbr addseq=random;\n")
	fcat("savetrees file=gk10.tre;\n")
	fcat("set autoclose=yes;\n\n")
	fcat("quit;\n\n\n")

    fcat("END;\n")
    close(zz) #Close the connection

}

# This Function writes the simulated table to a NEXUS formatted file
# For use with PAUP but all uninformative (i.e. invar/autapo) are removed
# Filename is the name of the output file to be created
# Table is the alignment created by PhyloSim
# tree is the tree for which simulations are conducted on, we will write it in the PAUP block
write.nexus.PAUP.nouninformative <- function(filename, table, tree){
    
    tree.to.write<-write.tree(tree) #get the tree as a newick formatted object
    
    
    #Open a connection and write the file
    #This is based on the write.nexus.data.R function in "APE"
    zz<-file(filename,"w")
    "fcat" <- function (..., file = zz) #Function to make concatenation easier
    {
        cat(..., file = filename, sep = "", append = TRUE)
    }
    indent<- "  "
    ntaxa<-length(table[,1])
    nchar<-length(table[1,])
    fcat("#NEXUS\n\n\n")
    NCHAR<- paste("NCHAR=", nchar, sep="")
    NTAX <- paste("NTAX=", ntaxa, sep = "")
    fcat("BEGIN DATA;\n")
    fcat("DIMENSIONS  ",NTAX,indent,NCHAR,";\n")
    fcat("FORMAT DATATYPE=STANDARD;\n")
    fcat("MATRIX\n\n\n\n\n\n")

    
    for (i in 1:length(table[,1])){    	
        fcat("\n", row.names(table)[i], indent, table[i,], "\n\n") 
    }
    
    
    fcat(indent, ";\n\n\n\n\n")

    fcat("END;\n")

    fcat("begin trees;\n\n\n")
    fcat("tree tree_1 = [&U]  ",tree.to.write,"\n\n")
    fcat(";\n\n")
    fcat("END;\n")


    fcat("begin paup;\n\n\n")

    fcat("exclude uninf;\n\n\n")

fcat("set increase=auto autoinc=100;\n")
fcat("hs nreps=1000 swap=tbr addseq=random;\n")
fcat("savetrees file=equal_Nouninf.tre;\n\n")

fcat("PSet goloboff=yes gk=2;\n") 
fcat("hs nreps=1000 swap=tbr addseq=random;\n")
fcat("savetrees file=gk2_Nouninf.tre;\n\n")

fcat("PSet goloboff=yes gk=3;\n") 
fcat("hs nreps=1000 swap=tbr addseq=random;\n")
fcat("savetrees file=gk3_Nouninf.tre;\n\n")

fcat("PSet goloboff=yes gk=6;\n") 
fcat("hs nreps=1000 swap=tbr addseq=random;\n")
fcat("savetrees file=gk6_Nouninf.tre;\n\n")

fcat("PSet goloboff=yes gk=10;\n")
fcat("hs nreps=1000 swap=tbr addseq=random;\n")
fcat("savetrees file=gk10_Nouninf.tre;\n")
fcat("set autoclose=yes;\n\n")
fcat("quit;\n\n\n")

    fcat("END;\n")
    close(zz) #Close the connection

}

form.shellscript <- function(filename,scriptname,directory,iteration,sample){
    
    part.dir<-paste(directory,"/",filename,sep="")
    part.dir<-paste("export RUNDIR=","\"",part.dir,"/",filename,"_",sample,"_",iteration,"\"", "\n",sep="")
    full.dir<-paste(directory,"/",filename,sep="")
    full.dir<-paste("export RUNFLAGS=","\"",full.dir,"/",filename,"_",sample,"/",filename,"_",sample,"_",iteration,"_MB.NEX","\"", "\n",sep="")
    
    zz<-file(scriptname,"w")
    
    "fcat" <- function (..., file = zz) #Function to make concatenation easier
    {
        cat(..., file = scriptname, sep = "", append = TRUE)
    }
    
    fcat("#!/bin/bash\n")
    fcat("#PBS -l nodes=1:ppn=8\n")
    fcat("#PBS -q veryshort\n")
    fcat("# Set the working directory for the job --------------------------------- \n")
    fcat(part.dir )
    fcat("# Name of application --------------------------------------------------- \n")
    fcat("export APPLICATION=\"/usr/local2/apps/MrBayes-3.2.2-MPI/bin/mb\"\n")
    fcat("# Any input files/switches\n")
    fcat(full.dir)
    fcat("# Change into the working directory -------------------------------------\n")
    fcat("cd $RUNDIR\n")
    fcat("# Generate the list of nodes the code will run on -----------------------\n")
    fcat("cat $PBS_NODEFILE\n")
    fcat("export nodes=`cat $PBS_NODEFILE`\n")
    fcat("export nnodes=`cat $PBS_NODEFILE | wc -l`\n")
    fcat("export confile=inf.$PBS_JOBID.conf\n")
    fcat("for i in $nodes; do\n")
    fcat("\techo ${i} >>$confile\n")
    fcat("done\n\n")
    fcat("# Execute the code ------------------------------------------------------\n")
    fcat("mpirun -np $nnodes -machinefile $confile $APPLICATION  $RUNFLAGS\n")
    
    close(zz) #Close the connection
}

#Function to form a filename to be passed to write nexus
form.filename <- function(filename, replicates){
	
	full.filename<-paste(filename, "_", replicates,  sep="")
	return(full.filename)
		}

#Function to create a parent dir filled with individual dirs
directory.creation <- function(filename, replicates, samples){
    #Make the parent dir
    original.wd <- getwd()
    dir.create(file.path(original.wd, filename))
    
    #make the individual dirs
    setwd(file.path(original.wd, filename))
    #Make the qsub folder
    dir.create(file.path(getwd(),"Qsubs"))    
    for (i in 1:samples){
        
        dir.create(file.path(getwd(),paste(filename,"_",i,sep="")))
        setwd(file.path(getwd(),paste(filename,"_",i,sep="")))
        
        for(j in 1:replicates){
        
        dir.create(file.path(getwd(),paste(filename,"_",i,"_",j,sep="")))
        
        }
        setwd(file.path(original.wd,filename))
        
    }
    
    #Return to top 
    setwd(original.wd)
}





