---
layout: article
date: 2017-08-05 11:46
title: "gmail内のメールを別のアカウントに移行する"
tags: [google-app-script, gmail]
comments: true
published: true
---

自分のアカウントのgmailの容量が増えてきて，google driveの容量を圧迫していたのでサブ垢にメールを全部移行する．

**注意: 本記事のメール移行スクリプトの実行は自己責任でお願いします．本記事を読んで実行した結果について，当方は責任を一切負いません．**

### 手順

1. 移行したい方のgoogle accountにログインし，google driveのページを開く．[Google Drive](https://drive.google.com/drive/my-drive)  
2. 左にあるNewと書いてあるボタンをクリックし，→More→Google Apps Scriptをクリックする．
![app-script-start](/assets/images/gmail-data-migration/app-script-start.png)  
注: Google Apps Scriptが選べないときは，New→other→Connect more appsでGoogle Apps Scriptを検索してインストールする．  
3. [Google Sheets](https://docs.google.com/spreadsheets/u/0/)に移動し，新しいシートを作成する．  
4. C2のセルに自分のメールアドレス，C3のセルに移行先のメールアドレスを入力して保存する．  
5. Google Apps Scriptのページにタブなどで戻り，以下の関数を作成，保存する．(デフォルトの関数は削除する)  

    ```
function gmailMigration(){
    var sheet = SpreadsheetApp.openById("<SHEET-ID>").getActiveSheet();
    var offset = 0;
    var limit = 100;
    var search =  sheet.getRange("C2").getValue()
    var fw_addr = sheet.getRange("C3").getValue()
  
    var threads = GmailApp.search(search, offset, limit);
  
    for (var i in threads) {
      var thread = threads[i];
      var msgs = thread.getMessages();
  
      for(var m in msgs){
        var msg = msgs[m];
        var from = msg.getFrom();
        var subject = msg.getSubject();
        var body = msg.getPlainBody().replace(/\r/g, "").replace(/\n\n/g, "\n").replace(/\n\n/g, "\n");
        var body = "From:" + from + "\n" + "Subj:" + subject + "\n\n" + body;
  
        GmailApp.sendEmail(fw_addr, subject, body);
      }
      Utilities.sleep(1000);
    }
}
    ```  
なお，\<SHEET-ID\>の部分は作成したシートのURLの"https://docs.google.com/spreadsheets/d/hogehoge/edit#gid=0"のうちのhogehogeの部分である．  
6. 三角矢印の実行ボタンをクリックする．  
7. クリック後，Authorization requiredのポップアップが出るので良ければReview Permissionを押す．  
![authorization-required](/assets/images/gmail-data-migration/authorization-required.png)  
8. スクリプトを続行するためのアカウントを聞かれるので，自分が今有効にしているアカウントを選択する．  
9. このスクリプトを本当に信頼するか，という旨のことを聞かれるので，よければ以下の手順でスクリプト作成画面に戻って実行する．  

ポップアップのAdvancedをクリック→Go to <hogehoge> (unsafe)をクリック→テキストの入力ボックスに"continue"と入力してNEXTボタンをクリック

これで最新のメールが500件分移行先のメールアドレスに届けられるようになった．offset, limitの値を変更すれば最新から見て何件目から初めて何件取得するかを指定できる．これを更に回す処理を書いてもよいが，どうせG Suiteやってない自分は一日100件しか送れないようなのでそこの実装は止めた．(Google Apps free editionは100件/日，G Suiteは1,500件/日)

### 参考
- [【Google Apps Script】 Gmail を一括転送するスクリプト(ラベルでもなんでも)](http://oryou.hatenablog.com/entry/2017/01/09/191913)

↑ソースコードはほとんどここのリンクのコピペです

- [GMailで受信したメールを翻訳して転送](http://paqalex.blog133.fc2.com/blog-entry-29.html)
- [Quotas for Google Services](https://developers.google.com/apps-script/guides/services/quotas#flexible_quotas_early_access)
