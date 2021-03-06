function [time_gauss,mv_length] = envelope_rect( current_t, frame_rate, amplitude )


mv_length = ceil( current_t * frame_rate );
time_gauss = amplitude * ones( 1, mv_length );