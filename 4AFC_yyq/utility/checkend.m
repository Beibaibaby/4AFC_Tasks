function checkend
[~,~,kc]=KbCheck;
if kc(KbName('escape'))
    ListenChar(0);
    ShowCursor;
    sca;
    error('�û���ֹ����');
end