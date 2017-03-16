%% as it is
close all 
par.datprocess.Twin.PsRF = [0 25];
par.datprocess.Twin.PsRF_lo = [0 25];
par.datprocess.Twin.SpRF = [-25 2];
par.datprocess.Twin.SpRF_lo = [-25 2];

predata1 = b3_FORWARD_MODEL_legacy20161204( model1,Kbase,par,trudata,id,0); predata10=predata1;
	
[ predata1 ] = predat_process( predata1,'PsRF',par);
[ predata1 ] = predat_process( predata1,'SpRF',par);
[ predata1 ] = predat_process( predata1,'PsRF_lo',par);
[ predata1 ] = predat_process( predata1,'SpRF_lo',par);

plot_TRUvsPRE(predata10,predata1)
clone_figure(57,61);
set(61, 'pos',[ 15 25 1915 529])

%% with the new forward mod

predata2 = b3_FORWARD_MODEL( model1,Kbase,par,trudata,id,0); predata20=predata2;
	
[ predata2 ] = predat_process( predata2,'PsRF',par);
[ predata2 ] = predat_process( predata2,'SpRF',par);
[ predata2 ] = predat_process( predata2,'PsRF_lo',par);
[ predata2 ] = predat_process( predata2,'SpRF_lo',par);

plot_TRUvsPRE(predata20,predata2)
clone_figure(57,62);
set(62, 'pos',[ 15 25 1915 529])

%% original comparison
plot_TRUvsPRE(predata10,predata20)
clone_figure(57,56);
plot_TRUvsPRE(predata1,predata2)
