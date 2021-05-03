using DrWatson
@quickactivate

using AlertPushover, ProgressMeter, Alert

pushover_alert!(token = "atf8nck7kgbmasra3dwv17sotcmxy4", user = "ueskshjvnfy8xsmqnt2gjhkbfbszmy")
alert("Your julia script is finished!")
