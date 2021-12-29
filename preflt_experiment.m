%train a model on 1, 5, 7 preflt
% get data from the next week
% run inference of 1, 5, 7 on it
% plot inferences against actual
dataFolder = "C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\";
filePattern = fullfile(dataFolder, '*.csv'); % Change to whatever pattern you need.
filenames = dir(filePattern);
model_names = {};
models = {};
X_data = {};
Y_data = {};
data_names = {};
for k = 1 : 5
    filename = filenames(k).name;
    fullFilename = fullfile(filenames(k).folder, filename);
    fullFilename = convertCharsToStrings(fullFilename);
    fprintf(1, 'Now reading %s\n', fullFilename);
    filter_lens = {1, 5, 7};
    for i = 1:3
        filter_len = filter_lens{i};
        [model, X, Y] = hacking_2021_07_20_loadalldatafromfile4training(fullFilename, ...
        "C:\Users\Victor\Desktop\fiber-draw\alldatatrainfromonefile", ...
        filter_len);
        X_data{end+1} = X;
        Y_data{end+1} = Y;
        models{end+1} = model;
        model_names{end+1} = fullFilename + "_model_" + int2str(filter_len);
        data_names{end+1} = fullFilename + "_filtered_" + int2str(filter_len);
        close all;
    end
    close all;
end
error_matrix = zeros(15, 15);
for i=1:length(models)
    model = models{i};
    for j=1:length(X_data)
        inputs = X_data{j};
        inputs_batch = inputs';
        outputs = Y_data{j};
        predictions = {};
        mse_total = 0.0;
        for b=1:length(inputs_batch)
            model = resetState(model);
            [model, model_prediction] = predictAndUpdateState(model, inputs_batch(b));
            prediction = model_prediction{1}(1,:);
            target = outputs{b}(1,:); % selects bfd
            mse_total = mse_total + mse(prediction, target);
        end
        n_batches = length(inputs_batch)
        error_matrix(i,j) = mse_total / n_batches;
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