package org.gephi.spectralspatialisation;

import javax.swing.Icon;
import javax.swing.JPanel;
import org.gephi.layout.spi.Layout;
import org.gephi.layout.spi.LayoutBuilder;
import org.gephi.layout.spi.LayoutUI;
import org.openide.util.NbBundle;
import org.openide.util.lookup.ServiceProvider;

/**
 *
 * @author CERI
 */
@ServiceProvider(service = LayoutBuilder.class)
public class SpectralSpatialisationBuilder implements LayoutBuilder {

    private final SpectralSpatialisationLayoutUI ui = new SpectralSpatialisationLayoutUI();

    @Override
    public String getName() {
        return NbBundle.getMessage(SpectralSpatialisationBuilder.class, "name");
    }

    @Override
    public SpectralSpatialisation buildLayout() {
        return new SpectralSpatialisation(this);
    }

    @Override
    public LayoutUI getUI() {
        return ui;
    }

    private static class SpectralSpatialisationLayoutUI implements LayoutUI {

        @Override
        public String getDescription() {
            return NbBundle.getMessage(SpectralSpatialisation.class, "description");
        }

        @Override
        public Icon getIcon() {
            return null;
        }

        @Override
        public JPanel getSimplePanel(Layout layout) {
            return null;
        }

        @Override
        public int getQualityRank() {
            return 0;
        }

        @Override
        public int getSpeedRank() {
            return 0;
        }
    }
}
