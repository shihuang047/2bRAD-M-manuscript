Website for Network:
	#MetagenoNets: Inference and Insights for Microbial Association Networks
	https://web.rniapps.net/metagenonets

Data preprocessing:
	Parameter		Choice
	Normalization	Don't Normalize
	Prevalence	atleast 0.0001 of max
	Occurence	in atleast 10% samples
	Transformation	Don't Transform

Categorical networks:
      According Group of Society
	Parameters  	Propagate to modules
	Algorithm  	NAMAP w/ Spearman
	p-value  		0.05
	Iterations  	100
	Corr. cutoff	CRITICAL-R