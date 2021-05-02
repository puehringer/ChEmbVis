import * as React from "react";
import { InterpolationPage } from "./InterpolationPage";
import { Navbar, Nav, Form, Button } from "react-bootstrap";
import { EActiveTabs, ICollection, IInterpolatedParticle, IRegistry } from "./interfaces";
import { EmbeddingPage } from "./EmbeddingPage";
import { downloadFile } from "./utils";
import { getRegistry } from "./utils/api";
import { FileUploadModal } from "./components/FileUploadModal";

function App() {
  const [importFileModalShow, setImportFileModalShow] = React.useState<boolean>(false);
  const [activeTab, setActiveTab] = React.useState<EActiveTabs>(EActiveTabs.EMBEDDING);
  const [collections, _setCollections] = React.useState<ICollection[]>([]);
  const [registry, setRegistry] = React.useState<IRegistry | null>(null);
  const [interpolationStructures, setInterpolationStructures] = React.useState<string[]>([
    "NC1CC1C(=O)c1ccc2ccccc2c1",
    "O=C(CN1C(=O)CSc2ccc(S(=O)C3CC3)cc21)NCc1cccnc1"
  ]);

  React.useEffect(() => {
    getRegistry()
      .then((registry) => setRegistry(registry))
      .catch((e) => {
        console.error("Error initializing registry", e);
      });
  }, []);

  React.useEffect(() => {
  }, [collections]);

  const setCollections = React.useCallback((collections: ICollection[]) => {
    collections.forEach((collection) =>
      collection.data.forEach((p, i) => {
        p.index = i;
        p.collection = collection.name;
      })
    );

    _setCollections(collections);
  }, [_setCollections]);

  const interpolationCollection = React.useMemo(
    () => collections.find((c) => c.name === "Interpolated") as ICollection<IInterpolatedParticle> | undefined,
    [collections]
  );

  const setInterpolationCollection = React.useCallback(
    (collection: ICollection<IInterpolatedParticle>) => {
      _setCollections((collections) => [...collections.filter((c) => c.name !== "Interpolated"), collection]);
    },
    [_setCollections]
  );

  return (
    <div className="vh-100 d-flex flex-column">
      <Navbar collapseOnSelect expand="lg" bg="dark" variant="dark">
        <Navbar.Brand href="#home">CDDD Explorer</Navbar.Brand>
        <Navbar.Toggle aria-controls="responsive-navbar-nav" />
        <Navbar.Collapse id="responsive-navbar-nav">
          <Nav className="mr-auto">
            <Nav className="mr-auto">
              <Nav.Link
                href="#"
                active={activeTab === EActiveTabs.EMBEDDING}
                onSelect={() => setActiveTab(EActiveTabs.EMBEDDING)}
              >
                Embedding
              </Nav.Link>
              <Nav.Link
                href="#"
                active={activeTab === EActiveTabs.INTERPOLATION}
                onSelect={() => setActiveTab(EActiveTabs.INTERPOLATION)}
              >
                Interpolation
              </Nav.Link>
              <Nav.Link href="/datasets/" target="_blank" active={false}>
                Datasets
              </Nav.Link>
            </Nav>
          </Nav>
          <Form inline>
            <Button
              onClick={() => {
                setImportFileModalShow(true);
              }}
            >
              Import
            </Button>
            <Button
              className="ml-2"
              disabled={collections.length === 0}
              onClick={() => {
                downloadFile(collections, "export");
              }}
            >
              Export
            </Button>
          </Form>
        </Navbar.Collapse>
      </Navbar>
      <FileUploadModal
        open={importFileModalShow}
        setOpen={setImportFileModalShow}
        onSave={(value) => {
          if (value) {
            try {
              setCollections(JSON.parse(value));
              setImportFileModalShow(false);
            } catch (e) {
              console.error("Error parsing imported file.");
            }
          }
        }}
      />
      <div className="container-fluid mt-2" style={{ flex: 1, overflow: "auto" }}>
        <div className="row" style={{ height: "100%", overflow: "auto", position: "relative" }}>
          {activeTab === EActiveTabs.EMBEDDING ? (
            <EmbeddingPage
              registry={registry}
              collections={collections}
              setCollections={setCollections}
              interpolationStructures={interpolationStructures}
              setInterpolationStructures={setInterpolationStructures}
              setActiveTab={setActiveTab}
            />
          ) : null}
          {activeTab === EActiveTabs.INTERPOLATION ? (
            <InterpolationPage
              collection={interpolationCollection}
              setCollection={setInterpolationCollection}
              structures={interpolationStructures}
              setStructures={(structures: string[]) => setInterpolationStructures(structures)}
            />
          ) : null}
        </div>
      </div>
    </div>
  );
}

export default App;
