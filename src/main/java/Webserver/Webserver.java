package Webserver;

import org.eclipse.jetty.server.Server;

/**
 *  Webserver set and configure the server
 *
 * Created by antoine on 17/03/17.
 */
public class Webserver {

    public static void main (String[] args) throws Exception {
        Server server = new Server(8080);
        server.setHandler(new WebserverHandler());

        server.start();
        server.join();
    }
}
