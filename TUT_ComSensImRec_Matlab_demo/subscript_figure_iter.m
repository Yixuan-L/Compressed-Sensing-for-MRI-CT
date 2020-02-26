% subscript_figure_iter
if figure_every_n_seconds>=0
    if (now-time_counter(2))>(figure_every_n_seconds*1.157412771135569e-005)||k==1||k==k_final
        if k==1||k==k_final
            figure(figureHandle);
        end
        if gcf==figureHandle  %% figure
            subplot(2,3,1)
            imshow([theta_hat_0],[0 1]);
            if exist('theta','var')
                title(['$\hat{\theta}^{(0)}$, PSNR: ',num2str(err_theta_0,number_of_digits_fig)],'interpreter','latex','fontsize',figTitleFontSize)
            else
                title(['$\hat{\theta}^{(0)}$'],'interpreter','latex','fontsize',figTitleFontSize)
            end

            subplot(2,3,2)
            imshow([filtered_invT_y_hat_k],[0 1]);
            if exist('theta','var')
                title(['$\Phi(\mathcal{T}^{-1}(\hat{y}^{(k)}))$, PSNR: ',num2str(err_filtered(k),number_of_digits_fig)],'interpreter','latex','fontsize',figTitleFontSize)
            else
                title(['$\Phi(\mathcal{T}^{-1}(\hat{y}^{(k)}))$'],'interpreter','latex','fontsize',figTitleFontSize)
            end

            subplot(2,3,3)
            imshow([theta_hat_k],[0 1]);
            if exist('theta','var')
                title(['$\hat{\theta}^{(k)}$, PSNR: ',num2str(err_theta_hat(k),number_of_digits_fig)],'interpreter','latex','fontsize',figTitleFontSize)
            else
                title(['$\hat{\theta}^{(k)}$'],'interpreter','latex','fontsize',figTitleFontSize)
            end

            subplot(2,2,3)
            plot(std_excite(1:k),'r-'),grid on;
            set(gca,'Xlim',[0 k]);
            title(['std\{$\eta_{k}$\}$=$',num2str(std_excite(k),3)],'interpreter','latex','fontsize',figTitleFontSize)

            subplot(2,2,4)
            plot(err_excitation(1:k),'r--'),grid on;
            if exist('theta','var')
                hold on
                plot(err_theta_hat(1:k),'m'),grid on;
                plot(err_filtered(1:k),'b'),grid on;
                hold off
                legend_handle=legend('PSNR exc.noise','PSNR $\hat{\theta}^{(k)}$','PSNR $\Phi(\mathcal{T}^{-1}(\hat{y}^{(k)}))$');
            else
                legend_handle=legend('PSNR exc.noise');
            end
            %             set(legend_handle,'interpreter','latex','location','best','orientation','vertical') %% this caused warnings
            set(legend_handle,'interpreter','latex','location','SouthOutside','orientation','horizontal')
            set(gca,'Xlim',[0 k]);
            title([ 'Iter. $k=$', num2str(k),  '  ($k_{\rm{final}}=$',num2str(k_final),')'],'interpreter','latex','fontsize',figTitleFontSize)
            drawnow
            time_counter(2)=now;
        end
        pause(eps)
    end
end