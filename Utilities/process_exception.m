function err_msg = process_exception(ME)
    err_txt_start = 'Error: %s (Error ID: "%s")\n\n';
    err_txt_stack_component = ['Error in: %s\n', ...
                               '  line %s: %s'];
    
    txt = err_txt_start;
    msgelements = [string(ME.message), string(ME.identifier)];
    
    element_filenames = {ME.stack.file};
    element_names = {ME.stack.name};
    element_lines = [ME.stack.line];

    n_stackelements = numel(ME.stack);
    for current_element = 1:n_stackelements
        if current_element==1
            txt = [txt,err_txt_stack_component];
        else
            txt = [txt,'\n',err_txt_stack_component];
        end
        name = '';
        for inner_idx = n_stackelements:-1:current_element
            if inner_idx==n_stackelements
                name = [name, element_names{inner_idx}];
            else
                name = [name, '>', element_names{inner_idx}];
            end
        end
        S = readlines(element_filenames{current_element});
        codeline = S(element_lines(current_element));
        msgelements = [msgelements, string(name), element_lines(current_element), string(codeline)];
    end

    err_msg = sprintf(txt, msgelements);
end