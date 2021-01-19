#From Johnason Hsu

#Cas9 screen:
docker run -v ${PWD}/:/DATA -w /DATA pinellolab/crisprsurf SURF_deconvolution -f [INSERT FILE NAME] -pert cas9 -sim_n 1000 -range 7
 
#dCas9 screen:
docker run -v ${PWD}/:/DATA -w /DATA pinellolab/crisprsurf SURF_deconvolution -f [INSERT FILE NAME] -pert cas9 -sim_n 1000 -range 20
