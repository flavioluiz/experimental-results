function run_save_test(name)
    global Tfinal;
    global output input times  T AMP1 AMP2 control_input;
    set_param(gcs,'SimulationCommand','connect');
    set_param(gcs,'SimulationCommand','start');
    Tfinal
    pause(Tfinal+20);
    outp = output.signals.values(:,:);
    outputvalues = outp;

    inp = input.signals.values(:,:);
    inputvalues = inp;   
    save(name,'inputvalues','outputvalues','times', 'control_input')
end