---
layout: article
date: 2017-07-29 21:00
title: "docker startでエラーが出た"
tags: [docker, centos, selinux]
comments: true
published: true
---

ec2インスタンスのdocker上でredmineを運用し始めたのだけど，
`docker start redmine`で起動させようとしたら以下のエラーが出てしばらくはまった．

```sh
Error response from daemon: devmapper: Error mounting '/dev/mapper/docker-202:1-12597455-b72245ffd57cc14382b64c23cad9581654a63a08dd4ccc5f8d088623b0c5f201' on '/var/lib/docker/devicemapper/mnt/b72245ffd57cc14382b64c23cad9581654a63a08dd4ccc5f8d088623b0c5f201': invalid argument
Error: failed to start containers: redmine
```

適当にググって見ると，他の人も同様の現象が起こっていた．
[Error response from daemon: devmapper: Error mounting '/dev/mapper/docker-253:2-8652374-' on '/var/lib/docker/devicemapper/mnt/': invalid argument](https://github.com/moby/moby/issues/29622)

SELinuxをdisableにした後，rebootすると上のエラーが出ると書いてある．
確かに前回の起動時にSELinuxをdisableにしていたので，自分も当てはまっている．


### SELinuxの設定をenforcingにする

仕方がないのでSELinuxをもう一度enforcingに戻すことにする．

```sh
$ sudo vim /etc/selinux/config
```
`SELINUX=disabled`と書いてある部分のdisabledをenforcingに変更する．
エディタ起動するのが面倒な人は以下のコマンドで．

```sh
$ sudo sed -i 's/SELINUX=disabled/SELINUX=enforcing/g' /etc/selinux/config
```

SELinux変更後はOSの再起動が必要なので，`sudo reboot`で再起動をする．
再起動後，`getenforce`コマンドで結果がEnforcingになっていれば変更完了．

### SELinuxが有効の状態で，Docker内のコンテナからmountファイルにアクセスができるようにする

しかし，この状態だとSELinuxのセキュリティのためにホスト側のディレクトリをdockerコンテナ側にmountしたVolumeにコンテナ側からアクセスができない．

```sh
$ sudo docker exec -it redmine bash
# redmineコンテナ内に入る
$ ls data-share-dir  # Permission denied
```

[Host volume mount fails when using SElinux and a mountpoint as source #23989](https://github.com/moby/moby/issues/23989)
によると，以下のコマンドでアクセスできるようだ．

```sh
$ chcon -Rt svirt_sandbox_file_t /path/to/docker-redmine/share/
```
でSELinuxあった状態でもアクセスできた．

```sh
$ sudo docker exec -it redmine bash
# redmineコンテナ内に入る
$ ls data-share-dir  # 今度は成功！
```

`chcon`はSELinuxのセキュリティコンテクストを変更するコマンドらしい．
`-R`で指定したディレクトリを再帰的に変更，`-t`でセキュリティコンテキストのタイプを変更する．

セキュリティコンテキストが変更されたかどうかは以下のコマンドで確認する．

```sh
$ ls -Z
drwxrwxr-x. fhiyo fhiyo unconfined_u:object_r:svirt_sandbox_file_t:s0 share
```

左からmode, user, group, security context, file nameの順番．
