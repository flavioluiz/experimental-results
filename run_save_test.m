function run_save_test(name)
    global Tfinal;
    global output input times  T AMP1 AMP2 control_input T_STOP;
    name
    set_param(gcs,'SimulationCommand','connect');
    set_param(gcs,'SimulationCommand','start');
    Tfinal
    pause(T_STOP+5);
    set_param(gcs,'SimulationCommand','stop');
    pause(5);
    outp = output.signals.values(:,:);
    outputvalues = outp;

    inp = input.signals.values(:,:);
    inputvalues = inp;   
    save(name,'inputvalues','outputvalues','times', 'control_input')
    
end