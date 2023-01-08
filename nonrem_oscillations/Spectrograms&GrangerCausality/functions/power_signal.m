function P=power_signal(x)
% Computes power of signal
L=length(x);
P=(norm(x)^2)/L;
end