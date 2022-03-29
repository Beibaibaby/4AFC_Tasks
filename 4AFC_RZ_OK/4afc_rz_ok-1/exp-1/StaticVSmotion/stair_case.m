function [new_state, new_up_down] = stair_case( old_state, old_up_down, response, design, step_size )


% response 1 corrent  0 is incorrect.
% design

new_state = old_state;

if response == 0
    
    new_up_down(1) = old_up_down(1) + 1;
    new_up_down(2) = 0;
    
    if new_up_down(1) == design(1)
        new_state = old_state + step_size;
        new_up_down(1) = 0;
    end
    
    
else
    
    new_up_down(1) = 0;
    new_up_down(2) = old_up_down(2) + 1;
    
    if new_up_down(2) == design(2)
        
        new_state = old_state - step_size;
        new_up_down(2) = 0;
        
    end
    
end
