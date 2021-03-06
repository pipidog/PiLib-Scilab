// **** Purpose ****
// PiLab template generator
// **** variables ****
// [project_name]: 1x1, string
// <= the project name
// **** Version ****
// 06/05/2014 1st version
// 06/07/2014 erase 'hop-manual' function
// **** Comment ****
// This function can generate a template file, so you don't have to 
// remember all setting parameters

function PiLab_template(project_name,job_type)
    fid=mopen(project_name+'_'+string(job_type)+'.plb','w');
    del_ind=0;
    select job_type
    case 'lat'
        mputl('lat.Const=[]                      '...
        +'// lattice constant, 1x1 real',fid);
        mputl('lat.Primitive=[]                  '...
        +'// Primitive vectors, (3x3/2x2/1x1)',fid);
        mputl('lat.Sublatt=[]                    '...
        +'// sublattice position, (nx3/nx2/nx1/)',fid);
        mputl('lat.VecOrder=[2]                      '...
        +'// Order of searching primitive vectros, 1x1 integer',fid);
        mputl('lat.NNCriterion=[0.2]                      '...
        +'// criterion on bond difference for an order, 1x1 integer',fid);
        mputl('lat.NNOrder=[2]                      '...
        +'// Order of printed neighbor orders, 1x1 integer',fid);       
    case 'hop'
        mputl('hop.SiteOrb=[]                    '...
        +'// specify orbital of each site, nx2, [site, l]',fid);
        mputl('hop.Order=[]                      '...
        +'// order of nearest coupling, 1x1 int, <= lat.Order',fid);
        mputl('hop.SKint=[]                      '...
        +'// SK parameters, nx7, [Orb1,Orb2,nn_order,ts,tp,td,tf]',fid);
        mputl('hop.LS=[0]                        '...
        +'// strength of LS coupling, 1x1, real', fid);
        mputl('hop.Filter=[10^-2]                '...
        +'// filiter of small hopping elements, 1x1, real',fid);              
        mputl('hop.Basis=[''c'']                    '...
        +'//  Basis of matrix, ''c'', ''s'', ''rc'', ''rs''',fid);
        mputl('hop.SelState=[]                   '...
        +'// input state labels to select states, 1xn, integer',fid);
        mputl('hop.OnsiteE=[]                    '...
        +'// Onsite energy of selected states, 1xn, real', fid);
    case 'scc'
        mputl('scc.HubU=[]                       '...
        +'// U for each state, [state_label, U] or blank', fid) 
        mputl('scc.Charge=[]                     '...
        +'// charge of each state, 1x total state',fid);
        mputl('scc.Mixing=[1]                    '...
        +'// mixing parameter, 0~1',fid);
        mputl('scc.Iteration=[30]                '...
        +'// maximal iterations, integer',fid);
        mputl('scc.Converge=[10^-3]              '...
        +'// convergence criterion, real, at least < 0.1',fid);
        mputl('scc.Mesh=[]                       '...
        +'// k-space mesh for Ef, (1x1,1x2,1x3), large for metal',fid);
        mputl('scc.Temperature=[0.1]             '...
        +'// smearing temperature',fid);
        mputl('scc.Memory=[200];                   '...
        +'// size for each buffer file(MB), suggest: 5~200',fid);
    case 'ban'
        mputl('ban.Format=[''coefficient'']      '...
        +'  // ''coefficient'' or ''coordinate''',fid);
        mputl('ban.Path=[]                       '...
        +'// points to defined your paths, nx3/nx2/nx1',fid);
        mputl('ban.DivType=''unit''              '...
        +'// how to divide k_path, ''unit'' or ''all''', fid);
        mputl('ban.Div=[20]                      '...
        +'// k points mesh of each path',fid);
        mputl('ban.Conv=''full''                 '...
        +'  // method to perform FT, ''full'' or ''cell''',fid);
        mputl('ban.Draw=[''on'']                 '...
        +'  // draw band structure, ''on'' or ''off''',fid); 
        mputl('ban.Shift=[''on'']                '...
        +'  // shift Ef to 0 in band plot, ''on'' or ''off',fid)
        mputl('ban.Ebound=[]                     '...
        +'  // energy bound of band plot',fid);
    case 'ban_den'
        mputl('ban_den.State=[]                  '...
        +'// specify band states, [k_label, eigenstate_label]',fid);
        mputl('ban_den.Draw=[''on'']               '...
        +'// whether to draw results',fid);
        mputl('ban_den.DrawType=[''site'']         '...
        +'// site/orb/suborb/site_up/site_dn/orb_up/orb_dn',fid);
    case 'dsa'
        mputl('dsa.Mesh=[]                       '...
        +'// k-space mesh,1x3/1x2/1x1',fid);
        mputl('dsa.Interval=[0.1]                '...
        +'// energy interval, 1x1, real',fid);
        mputl('dsa.Draw=[''on'']                 '...
        +'  // whether draw DOS, ''on'' or ''off''',fid);
        mputl('dsa.Shift=[''on'']                '...
        +'  // whether shift Ef=0 in plot, ''on'' or ''off''',fid)
    case 'chn'
        mputl('chn.Mesh=[]                       '...
        +'// k-space mesh of half BZ, 1x2, int',fid);
        mputl('chn.OccBand=[]                    '...
        +'// number of occupied bands, 1x1 int',fid);
        mputl('chn.Kdiff=[1,2]                   '...
        +'// Differential vector to avoid divergence',fid);
        mputl('chn.TRSplit=3*1e-4                '...
        +'// strength of TR splitting',fid);
        mputl('chn.PRSplit=5*1e-4                '...
        +'// strength of PR splitting',fid);
        mputl('chn.CFSplit=7*1e-4                '...
        +'// strength of CF splitting',fid);
    case 'zti'
        mputl('zti.Mesh=[]                       '...
        +'// mesh of half 2D BZ plane, 1x2 int',fid);
        mputl('zti.OccBand=[]                    '...
        +'// # of bands below Ef, int',fid);
        mputl('zti.Draw=[''off'']                '...
        +'  // Whether to draw n-field, ''on'' / ''off''',fid);
        mputl('zti.CFSplit=[3*1e-4]              '...
        +'  // Strength of crystal field breaking',fid);
        mputl('zti.PRSplit=[5*1e-4]              '...
        +'  // strength of parity breaking',fid);
    case 'flq'
        mputl('flq.Frequency=[]                  '...
        +'// field frequency, 1x1 real',fid);
        mputl('flq.Order=[]                      '...
        +'// order of photon process, 1x1, int',fid);
        mputl('flq.Amplitude=[]                  '...
        +'// AC amplitude, 1x1 / 1x2 / 1x3, real',fid); 
        mputl('flq.Phase=[]                      '...
        +'// AC phase 1x1 / 1x2 / 1x3, real',fid);
        mputl('flq.Mesh=[]                       '...
        +'// BZ mesh for Ef',fid);
        mputl('flq.Charge=[]                     '...
        +'// Floquet electrons per unit cell',fid);
        mputl('flq.Temp=[0.01]                     '...
        +'// smearing temperature for Ef',fid);
    case 'flq_ban'
        mputl('flq_ban.Format=[''coefficient'']  '...
        +'  // ''coefficient'' or ''coordinate''',fid);
        mputl('flq_ban.Path=[]                   '...
        +'// points to defined your paths, nx3/nx2/nx1',fid);
        mputl('flq_ban.DivType=[''unit'']        '...
        +'  // how to divide k_path, ''unit'' or ''all''', fid);
        mputl('flq_ban.Div=[20]                  '...
        +'// k points of each path',fid);
        mputl('flq_ban.Conv=''full''             '...
        +'  // method to perform FT, ''full'' or ''cell''',fid);
        mputl('flq_ban.Draw=[''on'']             '...
        +'  // draw band structure, ''on'' or ''off''',fid); 
        mputl('flq_ban.Shift=[''on'']            '...
        +'  // shift Ef to 0 in band plot, ''on'' or ''off',fid)
        mputl('flq_ban.Ebound=[]                     '...
        +'  // energy bound of band plot',fid)
    case 'flq_ban_den'
        mputl('flq_ban_den.State=[]                  '...
        +'// specify band states, [k_label, eigenstate_label]',fid);
        mputl('flq_ban_den.Draw=[''on'']               '...
        +'// whether to draw results',fid);
        mputl('flq_ban_den.DrawType=[''site'']         '...
        +'// site/orb/suborb/site_up/site_dn/orb_up/orb_dn',fid);
    case 'flq_chn'
        mputl('flq_chn.Mesh=[]                   '...
        +'// mesh of half-BZ, 1x2 int',fid);
        mputl('flq_chn.OccBand=[]                '...
        +'// number of occupied band',fid); 
        mputl('flq_chn.Kdiff=[1,2]             '...
        +'// Differential vector to avoid divergence',fid);   
    else 
        disp('Error: PiLab_create, No such job type!');
        del_ind=1;
    end
    mclose(fid);
    if del_ind==1 then
        mdelete(project_name+'_'+string(job_type)+'.plb');
    end
endfunction
