import * as React from "react";
import { Modal } from "react-bootstrap";
import { Button } from "react-bootstrap";
import { IEnabledProjection } from "../interfaces";
import { ProjectionSettingsForm } from "./form/ProjectionSettingsForm";

export function ProjectionSettingsModal({
  config: oldConfig,
  availableProjections,
  availableNearestNeighbors,
  availableClusters,
  availableProperties,
  availableOpacityProperties,
  availableConnectByProperties,
  onHide,
  onSave,
}: {
  config: IEnabledProjection | undefined;
  availableProjections?: string[];
  availableNearestNeighbors: string[];
  availableClusters: string[];
  availableProperties: string[];
  availableOpacityProperties: string[];
  availableConnectByProperties: string[];
  onHide?(): void;
  onSave?(value: IEnabledProjection): void;
}) {
  const [config, setConfig] = React.useState<IEnabledProjection | undefined>(oldConfig);

  React.useEffect(() => {
    setConfig(oldConfig);
  }, [oldConfig]);

  return (
    <>
      <Modal show={Boolean(config)} onHide={onHide} size="lg">
        <Modal.Header closeButton>
          <Modal.Title>Configure {config?.value}</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          {config ? (
            <>
              <div className="d-flex position-relative">
                <label
                  className={`my-1 me-2 col-sm-2 d-flex align-items-center`}
                  style={{ whiteSpace: "nowrap" }}
                  htmlFor="projectionNameInput"
                >
                  Name
                </label>
                <div className={`my-1 me-sm-2 flex-fill col-sm-10`}>
                  <input
                    type="text"
                    id="projectionNameInput"
                    className="form-control"
                    value={config.label}
                    onChange={(e) =>
                      setConfig({ ...config, value: e.currentTarget.value, label: e.currentTarget.value })
                    }
                  />
                </div>
              </div>
              <ProjectionSettingsForm
                inline={false}
                availableProjections={availableProjections}
                availableNearestNeighbors={availableNearestNeighbors}
                availableClusters={availableClusters}
                availableConnectByProperties={availableConnectByProperties}
                availableProperties={availableProperties}
                availableOpacityProperties={availableOpacityProperties}
                setConfig={setConfig}
                config={config}
              />
            </>
          ) : null}
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={onHide}>
            Close
          </Button>
          <Button variant="primary" onClick={() => (config ? onSave?.(config) : null)}>
            Save
          </Button>
        </Modal.Footer>
      </Modal>
    </>
  );
}
