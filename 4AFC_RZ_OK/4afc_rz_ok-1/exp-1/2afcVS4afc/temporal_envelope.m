function [profile,mv_length] = temporal_envelope(flat_duration, tail_duration, frame_rate, contrast)
% falt_duration + tail_duration is the total duration.
% durations are in sec. 
% 0.5*contrast is the peak of this temporal envelope.

flat_n_frames = round(frame_rate*flat_duration);
tail_n_frames = round(frame_rate*tail_duration);

if flat_n_frames > 0
    flat_profile = ones(1, flat_n_frames);
end

if tail_n_frames > 1 
    tail_frames = 1:tail_n_frames;
    temp_mean = mean(tail_frames);
    temp_sd = tail_n_frames/5;
    tail_profile = exp(-0.5*((tail_frames-temp_mean)/temp_sd).^2);
    head_part = tail_profile(1:floor(temp_mean));
    tail_part = tail_profile(floor(temp_mean)+1:end);
end


if flat_n_frames > 0
    if tail_n_frames > 1
        profile = [head_part, flat_profile, tail_part];
    else
        profile = flat_profile;
    end
else
    if tail_n_frames > 1
        profile = tail_profile;
    else
        profile = 0;
    end
end
        
profile = contrast*profile;
mv_length = length(profile);