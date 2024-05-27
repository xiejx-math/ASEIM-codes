function rhok=compute_rhok(errors)

k=length(errors)-1;
rhok=nthroot(errors(end),k);
end
