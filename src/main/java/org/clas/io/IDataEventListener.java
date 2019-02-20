package org.clas.io;

import org.jlab.io.base.DataEvent;

/**
 *
 * @author gavalian
 */
public interface IDataEventListener {
    void  dataEventAction(DataEvent event);
    void  timerUpdate();
    void  resetEventListener();
}
