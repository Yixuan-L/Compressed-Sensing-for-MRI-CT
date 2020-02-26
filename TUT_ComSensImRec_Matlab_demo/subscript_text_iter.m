%% subscript_text_iter.m

%% compute errors
err_excitation(k)=20*log10(1/(sqrt(1-numel(Omega)/numel(theta_hat_k))*std_excite(k)));   %%% PSNR IF THERE WERE ONLY NOISE
if exist('theta','var')    %% errors
    err_filtered(k)=10*log10(1/(mean((filtered_invT_y_hat_k(:)-theta(:)).^2)));          %%% PSNR
    err_theta_hat(k)=10*log10(1/(mean((theta_hat_k(:)-theta(:)).^2)));                   %%% PSNR
end

%% printout text
if text_every_n_seconds>=0
    if (now-time_counter(1))>(text_every_n_seconds*1.157412771135569e-005)||k==1||k==k_final
        if exist('theta','var')    %% errors
            disp(['Iter. k=',num2str(k),'  ','PSNR_filt: ',num2str(err_filtered(k),number_of_digits_text),'  std_exc: ',num2str(std_excite(max(1,k-1)),3),'  std_filt: ',num2str(std_filt,3),   '  ~',num2str(round((k_final/k-1)*(now-time_counter(3))/1.157412771135569e-005)),'s left']);
        else
            disp(['Iter. k=',num2str(k),'  ','  std_exc: ',num2str(std_excite(max(1,k-1)),3),'  std_filt: ',num2str(std_filt,3),   '  ~',num2str(round((k_final/k-1)*(now-time_counter(3))/1.157412771135569e-005)),'s left']);
        end
        time_counter(1)=now;
    end
end