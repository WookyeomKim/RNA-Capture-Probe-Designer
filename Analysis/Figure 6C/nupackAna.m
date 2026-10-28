
function nupackResult= nupackAna(material, strandA,Acon,strandB,Bcon, temp, sodiumC)

% Import the NUPACK module
nupack = py.importlib.import_module('nupack');

% Define physical model
my_model = nupack.Model(pyargs('material', material, 'celsius', temp, 'sodium', sodiumC));

% Define strand species
a = nupack.Strand(strandA, pyargs('name', 'a'));
% b = nupack.Strand('CAGTTCCTTGTCTATCGGC', pyargs('name', 'b'));

b= nupack.Strand(strandB,pyargs('name','b'));

% Create a Python dictionary for strands with their concentrations
strands = py.dict({{a, py.float(Acon)}, {b, py.float(Bcon)}});

% Define complexes
complexes = nupack.SetSpec(pyargs('max_size', int32(2)));

% Define tube ensemble containing strands at specified concentrations
t1 = nupack.Tube(pyargs('strands', strands, 'complexes', complexes, 'name', 'Tube t1'));

% Analyze the tube ensemble
compute = py.list({'pairs', 'mfe'});
tubes = py.list({t1});
tube_result = nupack.tube_analysis(pyargs('tubes', tubes, 'compute', compute, 'model', my_model));

% Access the 'complexes' attribute
complex_results = py.getattr(tube_result, 'complexes');

% Display available keys in the complexes dictionary
complex_keys = py.list(complex_results.keys());


% Extract delta G (dG) and MFE values from each complex
deltaG_values = [];
MFE_values=[];
for i = 1:length(complex_keys)
    key = complex_keys{i};
    complex_info = complex_results{key};
    deltaG = complex_info.free_energy; % Assuming dG is an attribute of the complex_info
    
    deltaG_values = [deltaG_values; double(deltaG)];
    MFE_values=[MFE_values; complex_info.mfe_stack];
end

tube_results = py.getattr(tube_result, 'tubes');
tube_keys = py.list(tube_results.keys());
concentration = [];
for i = 1:length(tube_keys)
    tube_key = tube_keys{i};
    tube_info = tube_results{tube_key};
    concentrations = py.getattr(tube_info, 'complex_concentrations');
    conc_keys = py.list(concentrations.keys());
    
    for j = 1:length(conc_keys)
        conc_key = conc_keys{j};
        conc_value = concentrations{conc_key};
        concentration = [concentration; double(conc_value)];
    end
end

% Extract equilibrium concentration
tube_results = py.getattr(tube_result, 'tubes');
tube_keys = py.list(tube_results.keys());
concentration = [];
for i = 1:length(tube_keys)
    tube_key = tube_keys{i};
    tube_info = tube_results{tube_key};
    concentrations = py.getattr(tube_info, 'complex_concentrations');
    conc_keys = py.list(concentrations.keys());
    

    for j = 1:length(conc_keys)
        conc_key = conc_keys{j};
        conc_value = concentrations{conc_key};
        concentration = [concentration; double(conc_value)];
    end
end

%Result order 
%(b+b)
%(a+a)
%(a)
%(a+b)
%(b)

%but concentration is at 3rd from the top
%(a)
%(a+a)
%(a+b)
%(b)
%(b+b)
% Display the extracted delta G values
deltaG_values = [deltaG_values(3); deltaG_values(5);deltaG_values(2);deltaG_values(4);deltaG_values(1)];
MFE_values= [MFE_values(3); MFE_values(5); MFE_values(2); MFE_values(4); MFE_values(1)];
concentration=[concentration(1); concentration(4); concentration(2);concentration(3); concentration(5)];
nupackResult=[deltaG_values MFE_values concentration];


end
