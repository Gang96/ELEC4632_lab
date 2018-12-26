figure(1)
subplot(2,1,1)
plot(ScopeData(:,1), ScopeData(:,2),'r', ScopeData(:,1),ScopeData(:,3), 'b');
legend('Deadbeat Responce', 'Non-Deadbeat Responce');
grid on;
ylim([-1 1])