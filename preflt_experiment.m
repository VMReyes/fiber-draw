%train a model on 1, 5, 7 preflt
% get data from the next week
% run inference of 1, 5, 7 on it
% plot inferences against actual
[model_1, X_1, Y_1] = hacking_2021_07_20_loadalldatafromfile4training("C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\DrawData_Tower48_2020-12-01_to2020-12-08.csv", ...
                                                          "C:\Users\Victor\Desktop\fiber-draw\alldatatrainfromonefile", ...
                                                          1);
save("alldatatrainfromonefile\model.mat", "model_1")
[model_5, X_5, Y_5] = hacking_2021_07_20_loadalldatafromfile4training("C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\DrawData_Tower48_2020-12-01_to2020-12-08.csv", ...
                                                          "C:\Users\Victor\Desktop\fiber-draw\alldatatrainfromonefile", ...
                                                          5);
[model_7, X_7, Y_7] = hacking_2021_07_20_loadalldatafromfile4training("C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\DrawData_Tower48_2020-12-01_to2020-12-08.csv", ...
                                                          "C:\Users\Victor\Desktop\fiber-draw\alldatatrainfromonefile", ...
                                                          7);
[model_test_1, XTest_1, YTest_1] = hacking_2021_07_20_loadalldatafromfile4training("C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\DrawData_Tower48_2020-12-08_to2020-12-15.csv", ...
                                                          "C:\Users\Victor\Desktop\fiber-draw\alldatatrainfromonefile", ...
                                                          1);
close all;
models = {model_1, model_5, model_7, model_test_1};
input_set = {X_1, X_5, X_7, XTest_1};
output_set = {Y_1, Y_5, Y_7, YTest_1};
error_matrix = zeros(4,4);
for i=1:length(models)
    model = models{i};
    for j=1:length(input_set)
        inputs = input_set{j};
        inputs_batch = inputs';
        outputs = output_set{j};
        
        model_prediction = model.predict(inputs_batch);
        
        model_prediction_transpose = model_prediction';
        
        error_matrix(i,j) = mse(model_prediction_transpose(:,1), outputs(:,1));
    end
end


% close all;
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'defaultTextInterpreter','latex');
% 
% y_pred = trainedNetworkModel.predict(XTest_1{1}');
% y_pred = y_pred';
% figure; plot(Y_sub(:,1)); hold on; plot(y_pred(:,1),'r'); 
% legend({'Data','Prediction'});
% 
% stl_plottrainingresults_FFTPSD_function(Y_sub(:,1), Y_sub(:,1), y_pred(:,1), y_pred(:,1));
% 
% figure(3000); xlim([0 1]); ylim([-5 4]); 
% figure(3001); grid minor;  